import numpy as np
import yaml

from classes import *
import functions as fc
import plots as pl

with open(r'configuration.yaml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    params = yaml.load(file, Loader=yaml.FullLoader)

#uInf = params['Ma1']*cInf
#mu = params['rhoInf']*uInf*params['C']/params['Re1']
#gcoeff = (gamma - 1) / 2.0 
#gpow = gamma / (gamma - 1)

#### Exit ###
objExit = DataSlice(params, 'exit')
ye = objExit.GetY()
ue = objExit.GetU()
ve = objExit.GetV()
pe = objExit.GetP()
rhoe = objExit.GetRho()
# Exit pressure
for i in range(params['nfiles']):
    P2 = fc.mass_flow_average_quantity(ye[i], ue[i], pe[i])
print('P2=',P2)

# Exit angle
for i in range(params['nfiles']):
    ubar = fc.mixed_out_average_quantity(ye[i], ue[i], params['Py']) # streamwise velocity
    vbar = fc.mixed_out_average_quantity(ye[i], ve[i], params['Py']) # normal velocity
alphaM = np.arctan(vbar/ubar)
uM = np.sqrt(ubar*ubar+vbar*vbar)
print('alphaM=', alphaM*180/np.pi)

# Exit Reynolds number
for i in range(params['nfiles']):
    rhobar = fc.mixed_out_average_quantity(ye[i], rhoe[i], params['Py']) # density
Re2 = (rhobar*params['C']*uM)/params['mu']
print('Re2=', Re2)

# Exit Mach number
for i in range(params['nfiles']):
    pbar = fc.mixed_out_average_quantity(ye[i], pe[i], params['Py']) # pressure
Tbar = pbar/(params['R']*rhobar)
cbar = np.sqrt(params['gamma']*params['R']*Tbar)
Ma2 = uM/cbar
print('Ma2=', Ma2)

### Inlet ###
objInlet = DataSlice(params, 'inlet')
yi = objInlet.GetY()
ui = objInlet.GetU()
rhoi = objInlet.GetRho()
vi = objInlet.GetV()
pi = objInlet.GetP()
for i in range(params['nfiles']):
    P1_s = fc.mass_flow_average_quantity(yi[i], ui[i], pi[i]) # Static pressure
    rhoA = fc.mass_flow_average_quantity(yi[i], ui[i], rhoi[i]) # density
    uA = fc.mass_flow_average_quantity(yi[i], ui[i], ui[i]) # streamwise velocity
    vA = fc.mass_flow_average_quantity(yi[i], ui[i], vi[i]) # normal velocity
P1 = P1_s + 0.5*rhoA*(uA*uA + vA*vA) 
print('P1=',P1)

# Pressure ratio
print('P2/P1=',P2/P1)

### Blade ###
# Initalisation: Data from .csv (wall)
objField = Data(params) # Class object
objExp = DataExp(params)

# Cp distribution
p = objField.GetP()
npoints = objField.GetNpoints()
cp = fc.cp_distribution(p, P1, P2, npoints, params)

# Axial chord
x = objField.GetX()
cax = fc.axial_chord(x[0]) # first file only

# Get experimental data
xexp = objExp.GetX()
cpexp = objExp.GetCP()

# Plot Cp distribution
pl.plot_cp(x/cax, cp, xexp, cpexp) 