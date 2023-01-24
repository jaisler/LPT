import numpy as np
import yaml

from classes import *
import functions as fc
import plots as pl

with open(r'configuration.yaml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    params = yaml.load(file, Loader=yaml.FullLoader)

#### Exit ###
objExit = DataSlice(params, 'exit')
ye = objExit.GetY()
ue = objExit.GetU()
ve = objExit.GetV()
pe = objExit.GetP()
rhoe = objExit.GetRho()
# Exit pressure
P2 = []
for i in range(params['nfiles']):
    P2.append(fc.mass_flow_average_quantity(ye[i], ue[i], pe[i]))

# Exit angle
alphaM = []
uM = []
for i in range(params['nfiles']):
    ubar = fc.mixed_out_average_quantity(ye[i], ue[i], params['Py']) # streamwise velocity
    vbar = fc.mixed_out_average_quantity(ye[i], ve[i], params['Py']) # normal velocity
    alphaM.append(np.arctan(vbar/ubar))
    uM.append(np.sqrt(ubar*ubar+vbar*vbar))

# Exit Reynolds number
Re2 = []
rhobar = []
for i in range(params['nfiles']):
    rhobar.append(fc.mixed_out_average_quantity(ye[i], rhoe[i], params['Py'])) # density
    Re2.append((rhobar[i]*params['C']*uM[i])/params['mu'][i])

# Exit Mach number
Ma2 = []
for i in range(params['nfiles']):
    pbar = fc.mixed_out_average_quantity(ye[i], pe[i], params['Py']) # pressure
    Tbar = pbar/(params['R']*rhobar[i])
    cbar = np.sqrt(params['gamma']*params['R']*Tbar)
    Ma2.append(uM/cbar)

### Inlet ###
objInlet = DataSlice(params, 'inlet')
yi = objInlet.GetY()
ui = objInlet.GetU()
rhoi = objInlet.GetRho()
vi = objInlet.GetV()
pi = objInlet.GetP()
P1_s = []
P1 = []
for i in range(params['nfiles']):
    P1_s.append(fc.mass_flow_average_quantity(yi[i], ui[i], pi[i])) # Static pressure
    rhoA = fc.mass_flow_average_quantity(yi[i], ui[i], rhoi[i]) # density
    uA = fc.mass_flow_average_quantity(yi[i], ui[i], ui[i]) # streamwise velocity
    vA = fc.mass_flow_average_quantity(yi[i], ui[i], vi[i]) # normal velocity
    P1.append(P1_s[i] + 0.5*rhoA*(uA*uA + vA*vA))

# Inlet Isentropic Mach Number (Mass average)
# Exit Isentropic Mach Number (Mass average)
gcoeff = 2.0 / (params['gamma'] - 1) 
gpow = (params['gamma'] - 1) / params['gamma']
Ma1is = []
Ma2is = []
for i in range(params['nfiles']):
    Ma1is.append(np.sqrt(((P1[i]/P1_s[i])**(gpow) - 1) * gcoeff))
    Ma2is.append(np.sqrt(((P1[i]/P2[i])**(gpow) - 1) * gcoeff))
 
 # Print data
for i in range(params['nfiles']):
    print('File'+str(i)+':')
    print('P2/P1=',P2[i]/P1[i])
    print('alphaM=', alphaM[i]*180/np.pi)
    print('Ma1is=', Ma1is[i])
    print('Ma2is=', Ma2is[i])
    print('Re2=', Re2[i])
    print('\n')

### Blade ### 
# Initalisation: Data from .csv (wall)
objField = Data(params) # Class object
objExp = DataExp(params) # Experimental data

# Cp distribution
p = objField.GetP()
npoints = objField.GetNpoints()
cp = fc.pressure_coefficient(p, P1, P2, npoints, params)

# Axial chord
x = objField.GetX()
cax = fc.axial_chord(x[0]) # first file only

# Get experimental data
xexp = objExp.GetX()
cpexp = objExp.GetCP()

# Plot Cp distribution
pl.plot_cp(x/cax, cp, xexp, cpexp, params['path0']) 

# X-shear stress 
wss = objField.GetWSS()
# Skin friction: Implementation from Garai et al. 2015
cf = fc.skin_friction_coefficiet(wss, P1, P2, params) 

# Plot x-shear stress
pl.plot_cf(x/cax, cf, params['path0'])