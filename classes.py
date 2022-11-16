import numpy as np
import pandas as pd
import functions as fc

class Data:
    def __init__(self, params):
        """ Initialisation of the data from .csv files. """
        df = []
        for i in range(params['nfiles']):
            df.append(pd.read_csv(params['path'+str(i)] + '/' 
                + params['fileF'] + '.csv', delimiter=','))
            
        self.x = [] 
        self.wss = []
        self.p = []
        self.npoints = []
        for i in range(params['nfiles']):
            self.x.append(df[i]['x'])
            self.p.append(df[i]['p'])
            self.wss.append(df[i]['Shear_s'])
            self.npoints.append(len(df[i]['x']))

    def GetX(self):       
        return self.x

    def GetP(self):
        return self.p

    def GetWSS(self):
        return self.wss

    def GetNpoints(self):
        return self.npoints

    def WritePressureToCSV(self, fileNumber, fileHandle):
        fileHandle.write('x, p' + '\n')
        for i in range(self.npoints[i]):
            fileHandle.write(str(self.x[fileNumber][i]) + ',' + str(self.p[fileNumber][i]) + '\n')

    def WriteShearStressToCSV(self, fileNumber, fileHandle):
        fileHandle.write('x, Shear_s' + '\n')
        for i in range(self.npoints):
            fileHandle.write(str(self.x[fileNumber][i]) + ',' + str(self.wss[fileNumber][i]) + '\n')

class DataInlet:
    def __init__(self, params):
        """ Initialisation of the data from .csv files. 
            The data in the file is related to the inlet
            measurements. It was taken at -0.5Cax from
            the leading edge """

        df = []
        for i in range(params['nfiles']):
            df.append(pd.read_csv(params['path0'] + '/' 
                + params['fileI'] + '.csv', delimiter=','))
            
        self.y = [] 
        self.p = []
        self.rho = []
        self.u = []
        self.v = []
        self.npoints = []
        for i in range(params['nfiles']):
            self.x.append(df[i]['x'])
            self.cp.append(df[i]['cp'])
            self.npoints.append(len(df[i]['x']))

    def GetNpoints(self):
        return self.npoints

    def GetX(self):       
        return self.x

class DataSlice:
    def __init__(self, params, loc):
        """ Initialisation of the data from .csv files. 
            The data in the file is related to the inlet or exit
            measurements. It is taken at -0.5*Cax or at Cax + 0.5Cax from
            the leading edge for the inlet and exit measurements, respectively 
            """
        
        df = []
        if loc == 'exit': 
            for i in range(params['nfiles']):
                df.append(pd.read_csv(params['path0'] + '/' 
                    + params['fileO'] + '.csv', delimiter=','))
        else:
            for i in range(params['nfiles']):
                df.append(pd.read_csv(params['path0'] + '/' 
                    + params['fileI'] + '.csv', delimiter=','))

        self.y = [] 
        self.ma = []
        self.rho = []
        self.p = []
        self.u = []
        self.v = []
        self.npoints = []
        for i in range(params['nfiles']):
            self.y.append(df[i]['y'])
            self.rho.append(df[i]['rho'])
            self.p.append(df[i]['p'])
            self.u.append(df[i]['u'])
            self.v.append(df[i]['v'])
            self.npoints.append(len(df[i]['y']))

        fc.bubble_sort(self.y, self.rho, self.p, self.u, self.v)

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

class DataExp:
    def __init__(self, params):
        """ Initialisation of the data from .csv files. """
        df = []
        for i in range(params['nfiles']):
            df.append(pd.read_csv(params['pathExp'] + '/' 
                + params['fileE'] + '.csv', delimiter=','))
            
        self.x = [] 
        self.cp = []
        self.npoints = []
        for i in range(params['nfiles']):
            self.x.append(df[i]['x'])
            self.cp.append(df[i]['cp'])
            self.npoints.append(len(df[i]['x']))
            
    def GetX(self):       
        return self.x

    def GetCP(self):
        return self.cp

    def GetNpoints(self):
        return self.npoints

