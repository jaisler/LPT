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

class DataTKE:
    def __init__(self, params):
        """ Extract data from a .csv file 
            to calculate the Power Density Spectra """

        self.y = []
        self.rho = []
        self.u = []
        self.v = []
        self.w = []
        self.p = []
        self.umag = []
        self.tau = []
        self.tke = []
        self.prod = []
        self.diss = []
        self.Tu = []
        self.t = []
        self.dt = []
        df = []
        for i in range(params['nfileTKE']):
            df.append(pd.read_csv(params['path'] + '/'
                + params['fileTKE'] + '.csv',
                delimiter=',',))

            # Prepare data: rho,u,v,w,p
            self.rho.append(df[i]['rho'])
            self.u.append(df[i]['rhou']/df[i]['rho'])
            self.v.append(df[i]['rhov']/df[i]['rho'])
            self.w.append(df[i]['rhow']/df[i]['rho'])
            self.umag.append(np.sqrt((df[i]['rhou']/df[i]['rho'])**2 
                                   + (df[i]['rhov']/df[i]['rho'])**2 
                                   + (df[i]['rhow']/df[i]['rho'])**2))
            self.p.append((params['gamma'] - 1) * (df[i]['E'] - 0.5 * self.umag[i]))

            # Calculate Reynolds average for rho
            rhom = np.zeros(params['npointsy'])
            l = 0
            tn = 0
            while(l < len(self.u[i])):
                for k in range(params['npointsy']):
                    rhoAvgZ = 0
                    for j in range(params['npointsz']):
                        rhoAvgZ += self.rho[i][j + k * params['npointsz'] + l]
                    # Obtain the average in z-direction
                    rhoAvgZ = rhoAvgZ / params['npointsz']
                    # Obtain the average in time
                    rhom[k] += rhoAvgZ 
                l += params['npointsz'] * params['npointsy']            
                tn += 1
            rhom /= tn 

            # Calculate Fabre averages for u,v,w -> {u_i} = <rho u_i> / <rho> 
            umf = np.zeros(params['npointsy'])
            vmf = np.zeros(params['npointsy'])
            wmf = np.zeros(params['npointsy'])
            l = 0
            tn = 0
            while(l < len(self.u[i])):
                for k in range(params['npointsy']):
                    uAvgZ = 0
                    vAvgZ = 0
                    wAvgZ = 0
                    for j in range(params['npointsz']):
                        uAvgZ += (self.rho[i][j + k * params['npointsz'] + l] 
                                  * self.u[i][j + k * params['npointsz'] + l]) 
                        vAvgZ += (self.rho[i][j + k * params['npointsz'] + l] 
                                  * self.v[i][j + k * params['npointsz'] + l])
                        wAvgZ += (self.rho[i][j + k * params['npointsz'] + l] 
                                  * self.w[i][j + k * params['npointsz'] + l])
                    # Obtain the average in z-direction
                    umf[k] += (uAvgZ / params['npointsz']) 
                    vmf[k] += (vAvgZ / params['npointsz']) 
                    wmf[k] += (wAvgZ / params['npointsz']) 
                l += params['npointsz'] * params['npointsy']            
                tn += 1
            umf /= (tn * rhom)
            vmf /= (tn * rhom)
            wmf /= (tn * rhom)

            # Calculate the Reynolds Stresses for compressible flow
            qunt = params['qunt']
            taumz = np.zeros((qunt, params['npointsy']))
            for j in range(params['npointsz']): # z
                tau = np.zeros((qunt, params['npointsy'])) # tau11: 0, tau22: 1, tau33: 2, tau12: 3
                tn = 0
                l = 0
                while(l < len(self.u[i])): # time
                    for k in range(params['npointsy']): # y
                        tau[0][k] += (self.rho[i][j + k * params['npointsz'] + l] 
                                   * (self.u[i][j + k * params['npointsz'] + l] - umf[k])**2)
                        tau[1][k] += (self.rho[i][j + k * params['npointsz'] + l]  
                                   * (self.v[i][j + k * params['npointsz'] + l] - vmf[k])**2)
                        tau[2][k] += (self.rho[i][j + k * params['npointsz'] + l]  
                                   * (self.w[i][j + k * params['npointsz'] + l] - wmf[k])**2)
                        tau[3][k] += (self.rho[i][j + k * params['npointsz'] + l] 
                                   * (self.u[i][j + k * params['npointsz'] + l] - umf[k]) 
                                   * (self.v[i][j + k * params['npointsz'] + l] - vmf[k]))
                    l += params['npointsy'] * params['npointsz']          
                    tn += 1

                for l in range(qunt):
                    tau[l] /= tn 
                    taumz[l] += tau[l]

            for i in range(qunt):
                taumz[i] /= params['npointsz'] 
                self.tau.append(taumz[i]) # Reynolds stresses

            # Calculate the TKE for compressible flow. 
            # tke is divided by rhom due to Favre average.
            self.tke.append((0.5 * (taumz[0] + taumz[1] + taumz[2])) / rhom) 

            # Turbulence Intensity units in percentage
            self.Tu.append(100 * (np.sqrt(0.333333*(taumz[0] + taumz[1] + taumz[2])) 
                           / params['Mach2is']))
            #self.prod.append(prod)
            #self.diss.append(diss)

            # Obtain y range based on the mesh
            #self.y = np.linspace(-1.378, -0.581, 80)
            self.y = np.linspace(-0.647, 0.151, 80)

    def GetTau(self):
        return self.tau

    def GetTKE(self):
        return self.tke 

    def GetProd(self):
        return self.prod   

    def GetDiss(self):
        return self.diss  

    def GetTu(self):
        return self.Tu  

    def GetY(self):
        return self.y

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
    
class DataPSDUpstream:
    def __init__(self, params):
        """ Extract data from a .csv file 
            to calculate the Power Density Spectra """

        self.u = [] 
        self.t = []
        self.dt = []
        df = []
        for i in range(params['nfilepsdUp']):
            df.append(pd.read_csv(params['path'] + '/'
                + params['filePSDUp'] + str(i) + '.csv',
                delimiter=',',
                skiprows=[j for j in range(1, params['psdlocUp'])]))

            # Obtain data in each psdNpoints rows
            df[i] = df[i].iloc[::params['psdNpointsUp'], :]

            self.u.append(df[i]['rhou']/df[i]['rho'])
            self.t.append(df[i]['t'])
            self.dt.append(df[i]['t'][params['psdNpointsUp']] - df[i]['t'][0])

    def GetU(self):
        return self.u

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

