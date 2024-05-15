import numpy as np
import pandas as pd
import functions as fc
from scipy.interpolate import CubicSpline
from scipy.interpolate import make_interp_spline

class Data:
    def __init__(self, params):
        """ Initialisation of the data from .csv files. """
        df = []
        self.x = [] 
        self.wss = []
        self.rho = []
        self.p = []
        self.npoints = []
        for i in range(params['nfiles']):
            df.append([])
            for j in range(params['nfileBlade']):
                df[i].append(pd.read_csv(params['path' + str(i)] + '/' 
                    + params['fileB'] + str(j) + '.csv', delimiter=','))

            wss = 0
            p = 0
            rho = 0
            for j in range(params['nfileBlade']):
                p = p + df[i][j]['p']
                wss = wss + df[i][j]['Shear_x']
                rho = rho + (df[i][j]['p'] / (params['R'] * df[i][j]['T'])) 

            self.x.append(df[i][j]['x'])
            self.p.append(p/params['nfileBlade'])
            self.wss.append(wss/params['nfileBlade'])
            self.rho.append(rho/params['nfileBlade'])
            self.npoints.append(len(df[i][j]['x']))

    def GetX(self):       
        return self.x

    def GetRho(self):
        return self.rho

    def GetP(self):
        return self.p

    def GetWSS(self):
        return self.wss

    def GetNpoints(self):
        return self.npoints
    
    def GetCf(self, params):
        cf = []
        for i in range(params['nfiles']):
            cf.append(self.wss[i]/(0.5*params['rhoInf']*params['uInf']**2))
        return cf

    def WritePressureToCSV(self, fileNumber, fileHandle):
        fileHandle.write('x, p' + '\n')
        for i in range(self.npoints[i]):
            fileHandle.write(str(self.x[fileNumber][i]) + ',' + str(self.p[fileNumber][i]) + '\n')

    def WriteShearStressToCSV(self, fileNumber, fileHandle):
        fileHandle.write('x, Shear_s' + '\n')
        for i in range(self.npoints):
            fileHandle.write(str(self.x[fileNumber][i]) + ',' + str(self.wss[fileNumber][i]) + '\n')

class DataSlice:
    def __init__(self, params, loc):
        """ Initialisation of the data from .csv files. 
            The data in the file is related to the inlet or exit
            measurements. It is taken at -0.3*Cax or at Cax + 0.5Cax from
            the leading edge for the inlet and exit measurements, respectively 
            """
        
        df = []
        if loc == 'exit': 
            for i in range(params['nfiles']):
                df.append([])
                for j in range(params['nfileExit']):
                    df[i].append(pd.read_csv(params['path' + str(i)] + '/' 
                        + params['fileO'] + str(j) + '.csv', delimiter=','))
                    nf = params['nfileExit']
        else:
            for i in range(params['nfiles']):
                df.append([])
                for j in range(params['nfileInlet']):
                    df[i].append(pd.read_csv(params['path' + str(i)] + '/' 
                        + params['fileI'] + str(j) + '.csv', delimiter=','))
                    nf = params['nfileInlet']

        # limits
        limit = []
        limit.append(df[0][0]['y'].min()) 
        limit.append(df[0][0]['y'].max()) 

        self.y = [] 
        self.rho = []
        self.p = []
        self.u = []
        self.v = []
        self.w = []
        self.npoints = []
        for i in range(params['nfiles']):
            y = 0
            rho = 0
            p = 0
            u = 0
            v = 0
            w = 0
            for j in range(nf): 
                fc.bubble_sort(df[i][j]['y'], df[i][j]['rho'], df[i][j]['p'], 
                               df[i][j]['u'], df[i][j]['v'], df[i][j]['w'])

                ysort = []
                rhosort = []
                usort = []
                vsort = []
                wsort = []
                psort = []
                # Clean vector: remove repeated values
                k = 0
                while k < (len(df[i][j]['y']) - 1):
                    ysort.append(df[i][j]['y'][k])
                    rhosort.append(df[i][j]['rho'][k])
                    usort.append(df[i][j]['u'][k])
                    vsort.append(df[i][j]['v'][k])
                    wsort.append(df[i][j]['w'][k])
                    psort.append(df[i][j]['p'][k])

                    if (df[i][j]['y'][k] >= df[i][j]['y'][k+1]):
                        k += 3
                    k += 1

                # Interpolation
                degree = 1
                csRho = make_interp_spline(ysort, rhosort, degree)
                csU = make_interp_spline(ysort, usort, degree)
                csV = make_interp_spline(ysort, vsort, degree)
                csW = make_interp_spline(ysort, wsort, degree)
                csP = make_interp_spline(ysort, psort, degree)

                # Generate 100 equispaced points along the x-axis from the minimum 
                # to the maximum x values
                yrange = np.linspace(limit[0], limit[1], params['spline'])

                # Points in a equispaced form
                rhoInt = csRho(yrange) 
                uInt = csU(yrange) 
                vInt = csV(yrange) 
                wInt = csW(yrange) 
                pInt = csP(yrange) 

                # Sum up for average
                y += yrange #+ df[i][j]['y']
                rho += rhoInt #+ df[i][j]['rho']
                p += pInt #+ df[i][j]['p']
                u += uInt #+ df[i][j]['u']
                v += vInt #+ df[i][j]['v']
                w += wInt #+ df[i][j]['w'] 
            
            # Spanwise average
            self.y.append(y/nf)
            self.rho.append(rho/nf)
            self.p.append(p/nf)
            self.u.append(u/nf)
            self.v.append(v/nf)
            if (params['spacedim'] == 3):
                self.w.append(w/nf)
            else:
                self.w.append(df[i][j]['w']*0.0) 
            self.npoints.append(len(df[i][j]['y']))

    def GetNpoints(self):
        return self.npoints

    def GetY(self):       
        return self.y

    def GetRho(self):       
        return self.rho

    def GetP(self):       
        return self.p

    def GetU(self):       
        return self.u

    def GetV(self):       
        return self.v

    def GetW(self):       
        return self.w

class DataPSD:
    def __init__(self, params):
        """ Extract data from a .csv file 
            to calculate the Power Density Spectra """

        self.u = [] # attribute or member
        self.v = []
        self.w = []
        self.umag = []
        self.p = []
        self.t = []
        self.dt = []
        df = []
        for i in range(params['nfilepsd']):
            #subprocess.call(["sed -e 's/ /,/g' "
            #          + params['path'] + "/"
            #          + params['filePSD'] + str(i)+ ".his > "
            #          + params['path'] + "/"
            #          + params['filePSD'] + str(i) + ".csv"], shell=True)
            df.append(pd.read_csv(params['path'] + '/'
                + params['filehist'] + str(i) + '.csv',
                delimiter=',',
                skiprows=[j for j in range(1, params['psdloc'])]))

            # Obtain data in each psdNpoints rows
            df[i] = df[i].iloc[::params['psdNpoints'], :]

            self.u.append(df[i]['rhou']/df[i]['rho'])
            self.v.append(df[i]['rhov']/df[i]['rho'])
            self.w.append(df[i]['rhow']/df[i]['rhou'])
            self.umag.append(np.sqrt(self.u[i]**2 + 
                self.v[i]**2 + self.w[i]**2))
            self.t.append(df[i]['t'])
            self.dt.append(df[i]['t'][params['psdNpoints']] - df[i]['t'][0])

    def GetU(self):
        return self.u

    def GetV(self):
        return self.v

    def GetW(self):
        return self.w

    def GetUmag(self):
        return self.umag

    def GetT(self):
        return self.t

    def GetDt(self):
        return self.dt

class DataSpatialCorrelation:
    def __init__(self, params):
        """ Extract data from a .csv file 
            to calculate the Power Density Spectra """

        self.u = [] # attribute or member
        self.v = []
        self.w = []
        self.t = []
        self.dt = []
        df = []
        for i in range(params['nfilecorr']):
            #subprocess.call(["sed -e 's/ /,/g' "
            #          + params['path'] + "/"
            #          + params['fileCorr'] + str(i)+ ".his > "
            #          + params['path'] + "/"
            #          + params['fileCorr'] + str(i) + ".csv"], shell=True)
            df.append(pd.read_csv(params['path'] + '/'
                + params['filehist'] + str(i) + '.csv',
                delimiter=',',))
            
            self.u.append(df[i]['rhou']/df[i]['rho'])
            self.v.append(df[i]['rhov']/df[i]['rho'])
            self.w.append(df[i]['rhow']/df[i]['rho'])
            self.t.append(df[i]['t'])
            self.dt.append(df[i]['t'][params['psdNpoints']]-df[i]['t'][0])

            # Mean velocity
            k = 0
            tn = 0
            denu = 0
            denv = 0
            denw = 0
            while (k < len(self.u[i])):
                meanuz = 0
                meanvz = 0
                meanwz = 0
                for j in range(params['npcorr']): 
                    meanuz += self.u[i][j+k]
                    meanvz += self.v[i][j+k]
                    meanwz += self.w[i][j+k]
                meanuz /= (params['npcorr'])
                meanvz /= (params['npcorr'])
                meanwz /= (params['npcorr'])
                denu += meanuz
                denv += meanvz
                denw += meanwz
                k += params['npcorr']
                tn += 1
            denu /= tn
            denv /= tn
            denw /= tn

            self.u[i] = self.u[i] - denu
            self.v[i] = self.v[i] - denv
            self.w[i] = self.w[i] - denw

    def GetUprime(self):
        return self.u

    def GetVprime(self):
        return self.v

    def GetWprime(self):
        return self.w

    def GetT(self):
        return self.t

    def GetDt(self):
        return self.dt

class DataExp:
    def __init__(self, params):
        """ Initialisation of the data from .csv files. """
        df = []
        for i in range(params['nfilesCp']):
            df.append(pd.read_csv(params['pathExp'] + '/' 
                + params['fileE'] + str(i) + '.csv', delimiter=','))
            
        self.x = [] 
        self.cp = []
        self.npoints = []
        for i in range(params['nfilesCp']):
            self.x.append(df[i]['x'])
            self.cp.append(df[i]['cp'])
            self.npoints.append(len(df[i]['x']))
            
    def GetX(self):       
        return self.x

    def GetCP(self):
        return self.cp

    def GetNpoints(self):
        return self.npoints

class DataExpCf:
    def __init__(self, params):
        """ Initialisation of the data from .csv files. """
        df = []
        for i in range(params['nfilesCf']):
            df.append(pd.read_csv(params['pathExpCf'] + '/' 
                + params['fileECf'] + str(i) + '.csv', delimiter=','))
            
        self.x = [] 
        self.cf = []
        self.npoints = []
        for i in range(params['nfilesCf']):
            self.x.append(df[i]['x'])
            self.cf.append(df[i]['cf'])
            self.npoints.append(len(df[i]['x']))
            
    def GetX(self):       
        return self.x

    def GetCF(self):
        return self.cf

    def GetNpoints(self):
        return self.npoints
    
class DataExpLoss:
    def __init__(self, params):
        """ Initialisation of the data from .csv files. """
        df = []
        for i in range(params['nfilesLoss']):
            df.append(pd.read_csv(params['pathExpLoss'] + '/' 
                + params['fileEL'] + str(i) +'.csv', delimiter=','))
            
        self.x = [] 
        self.loss = []
        self.npoints = []
        for i in range(params['nfilesLoss']):
            self.x.append(df[i]['x'])
            self.loss.append(df[i]['loss'])
            self.npoints.append(len(df[i]['x']))
            
    def GetX(self):       
        return self.x

    def GetLoss(self):
        return self.loss

    def GetNpoints(self):
        return self.npoints

