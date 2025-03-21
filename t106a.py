import numpy as np
import yaml

from classes import *
import functions as fc
import plots as pl

with open(r'configuration.yaml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    params = yaml.load(file, Loader=yaml.FullLoader)

if (params['routine']['bltu']['pre']):
    objTKEBL = DataTKEBL(params)

if (params['routine']['bltu']['post']):
    scale = 0.725
    BLtu = []
    tu = []
    up = []
    vp = []
    ut = []
    s = []
    for j in range(params['nfilesTKEBL']):
        BLtu.append([])
        tu.append([])
        up.append([])
        vp.append([])
        ut.append([])
        s.append([])
        for i in range(params['nfilesBLtu']):
            BLtu[j].append(fc.LoadCSV(params['path'+str(j)]+'/BLtupost'+str(i)+'.csv'))
            tu[j].append(BLtu[j][i]['tu'])
            up[j].append(BLtu[j][i]['up'])
            vp[j].append(BLtu[j][i]['vp'])
            s[j].append(BLtu[j][i]['s']* scale)
        ut[j] = (fc.mean_tangential_velocity(up[j], vp[j], params))
        fc.boundary_layer_adjustment(400, tu[j], params)
        fc.boundary_layer_adjustment(40, ut[j], params)

    pl.plot_tu_BL(s,tu,params['path0'])
    pl.plot_ut_BL(s,ut,params)

if (params['routine']['tke']['pre']):
    # Literature: Experiments and DNS
    objExpTau = DataExpTau(params)
    objExpTKE = DataExpTKE(params)
    yexptau = objExpTau.GetY()
    tauexp = objExpTau.GetTau()
    yexptke = objExpTKE.GetY()
    tkeexp = objExpTKE.GetTKE()

    # Nektar++
    objTKE = DataTKE(params)
    tau = objTKE.GetTau()
    tke = objTKE.GetTKE()
    Tu = objTKE.GetTu()
    y = objTKE.GetY()

    # Trubulence intensity Tu
    pl.plot_Tu(y,Tu,params['pathTKE'])

    # Reynolds Stresses
    ystar = fc.norm_position(y, params['yminTau'], params['Py'], params)
    (tauShift, ystarShift) = fc.shift_data(tau, ystar, params['shiftTau'])
    ystarFlip = fc.flip_data(ystarShift, params)
    pl.plot_tau(ystarFlip,tauShift,yexptau,tauexp,params['pathTKE'])
    
    # TKE
    yexptkestar = fc.norm_position(yexptke, params['yminTkeExp'], params['Py'], params)
    ystar = fc.norm_position(y, params['yminTau'], params['Py'], params)
    (tkeShift, ystarShift) = fc.shift_data(tke, ystar, params['shiftTke'])
    pl.plot_tke(ystarShift,tkeShift,yexptkestar,tkeexp,params['pathTKE'])

    del objTKE

if (params['routine']['tke']['post']):
    tu = []
    y = []
    df = []
    for i in range(params['nfilesTKEdns']):
        if (i == 0):
            add = -0.17
        if (i == 1):
            add = -0.08
        if (i == 2):
            add = -0.06
        
        df.append(fc.LoadCSV(params['path'+str(i)]+ '/' + params['fileTKE'] + '_post.csv'))
        tu.append(df[i]['tu'] + add)
        y.append(df[i]['s'] - 0.647)
            
    # Average:
    for i in range(params['nfilesTKEdns']):
        print("i = ", i, ", avg = ", np.average(tu[i]))

    pl.plot_Tu(y,tu,params['path0'])

if (params['routine']['psd']):
    objPSDUp = DataPSDUpstream(params)
    dt = objPSDUp.GetDt()
    u = objPSDUp.GetU()
    (fu, Su) = fc.CalculatePSD(u, dt, params, flag='up')
    for i in range(len(fu)):
        for j in range(len(fu[i])):
            fu[i][j] = fu[i][j] * (params['C'] / params['Uinf'])
            Su[i][j] = Su[i][j] / (params['Uinf']**2)
    pl.plot_psd_velocity_upstream(fu, Su, params)  
    # delete object  
    del objPSDUp

    objPSD = DataPSD(params)
    dt = objPSD.GetDt()
    u = objPSD.GetU()
    (fu, Su) = fc.CalculatePSD(u, dt, params, flag='down')
    for i in range(len(fu)):
        for j in range(len(fu[i])):
            fu[i][j] = fu[i][j] * (params['C'] / params['Uinf'])
            Su[i][j] = Su[i][j] / (params['Uinf']**2)
    pl.plot_psd_velocity(fu, Su, params)  
    # delete object  
    del objPSD

if (params['routine']['correlation']):
    objCorr = DataSpatialCorrelation(params)
    up = objCorr.GetUprime()
    vp = objCorr.GetVprime()
    wp = objCorr.GetWprime()

    Rxx = []
    Ryy = []
    Rzz = []
    for i in range(params['nfilescorr']):
        Rxx.append(fc.CalculateTwoPointsCorrelation(up[i][params['nfilecorr']-1], params))
        Ryy.append(fc.CalculateTwoPointsCorrelation(vp[i][params['nfilecorr']-1], params))
        Rzz.append(fc.CalculateTwoPointsCorrelation(wp[i][params['nfilecorr']-1], params))
    z = np.linspace(0, params['Lz'], num=params['npcorr'])
    pl.plot_tpcorr(z, Rxx, Ryy, Rzz, params['path0'])
    
    #Rii = []
    #Rii.append(fc.CalculateOnePointCorrelation(up[0], params))
    #Rii.append(fc.CalculateOnePointCorrelation(vp[0], params))
    #Rii.append(fc.CalculateOnePointCorrelation(wp[0], params))
    #t = np.linspace(0, params['time'], num=int(params['numSteps'] / 
    #                                        params['pfreq']) + 1) / params['tRef']
    #pl.plot_opcorr(t, Rii, params['path'])

    # delete object  
    del objCorr

if (params['routine']['phisical']):    
    #### Exit ###
    objExit = DataSlice(params, 'exit')
    ye = objExit.GetY()
    ue = objExit.GetU()
    ve = objExit.GetV()
    we = objExit.GetW()
    pe = objExit.GetP()
    rhoe = objExit.GetRho()
    # Exit pressure
    P2 = []
    for i in range(params['nfiles']):
        P2.append(fc.mass_flow_average_quantity(ye[i], ue[i], pe[i]))

    # Mixed-out exit angle
    alphaM = []
    uM = []
    for i in range(params['nfiles']):
        ubar = fc.mixed_out_average_quantity(ye[i], ue[i], params['Py']) # streamwise velocity
        vbar = fc.mixed_out_average_quantity(ye[i], ve[i], params['Py']) # normal velocity
        wbar = fc.mixed_out_average_quantity(ye[i], we[i], params['Py']) # normal velocity
        alphaM.append(np.arctan(vbar/ubar))
        uM.append(np.sqrt(ubar*ubar+vbar*vbar+wbar*wbar))

    # Mass average exit angle
    alphaA = []
    for i in range(params['nfiles']):
        ubar = fc.mass_flow_average_quantity(ye[i], ue[i], ue[i]) # streamwise velocity
        vbar = fc.mass_flow_average_quantity(ye[i], ue[i], ve[i]) # normal velocity
        alphaA.append(np.arctan(vbar/ubar))

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

    # Exit Reynolds number in each point
    Re2p = []
    for i in range(len(ye[0])):
        Re2p.append(rhoe[0][i]*params['C']*np.sqrt(ue[0][i]*ue[0][i]+ve[0][i]*ve[0][i]+
                    we[0][i]*we[0][i]) / params['mu'][0])   

    # Plot Reynolds   
    #pl.plot_Re(ye, Re2p, params['path'])  

    # Exit Ma number in each point
    Mach2p = []
    for i in range(len(ye[0])):
        Tp = pe[0][i] / (params['R'] * rhoe[0][i])
        Mach2p.append(np.sqrt(ue[0][i]*ue[0][i]+ve[0][i]*ve[0][i]+
                              we[0][i]*we[0][i]) / np.sqrt(params['gamma']*params['R'] * Tp))   

    # Plot Reynolds   
    #pl.plot_Mach(ye, Mach2p, params['path'])  

    ### Inlet ###
    objInlet = DataSlice(params, 'inlet')
    yi = objInlet.GetY()
    rhoi = objInlet.GetRho()
    ui = objInlet.GetU()
    vi = objInlet.GetV()
    wi = objInlet.GetW()
    pi = objInlet.GetP()
    P1_s = []
    P1 = []
    for i in range(params['nfiles']):
        P1_s.append(fc.mass_flow_average_quantity(yi[i], ui[i], pi[i])) # Static pressure
        rhoA = fc.mass_flow_average_quantity(yi[i], ui[i], rhoi[i]) # density
        uA = fc.mass_flow_average_quantity(yi[i], ui[i], ui[i]) # streamwise velocity
        vA = fc.mass_flow_average_quantity(yi[i], ui[i], vi[i]) # normal velocity
        wA = fc.mass_flow_average_quantity(yi[i], ui[i], wi[i]) # spanwise velocity

        P1.append(P1_s[i] + 0.5*rhoA*(uA*uA + vA*vA + wA*wA))

    # Inlet angle
    alphaMI = []
    uMI = []
    for i in range(params['nfiles']):
        ubar = fc.mixed_out_average_quantity(yi[i], ui[i], params['Py']) # streamwise velocity
        vbar = fc.mixed_out_average_quantity(yi[i], vi[i], params['Py']) # normal velocity
        wbar = fc.mixed_out_average_quantity(yi[i], wi[i], params['Py']) # normal velocity
        alphaMI.append(np.arctan(vbar/ubar))
        uMI.append(np.sqrt(ubar*ubar+vbar*vbar+wbar*wbar))

    # Inlet Isentropic Mach Number (Mass average)
    # Exit Isentropic Mach Number (Mass average)
    gcoeff = 2.0 / (params['gamma'] - 1) 
    gpow = (params['gamma'] - 1) / params['gamma']
    Ma1is = []
    Ma2is = []
    for i in range(params['nfiles']):
        Ma1is.append(np.sqrt(((P1[i]/P1_s[i])**(gpow) - 1) * gcoeff))
        Ma2is.append(np.sqrt(((P1[i]/P2[i])**(gpow) - 1) * gcoeff))

    # Exit Ma number in each point
    Mach1p = []
    for i in range(len(yi[0])):
        Tp = pi[0][i] / (params['R'] * rhoi[0][i])
        Mach1p.append(np.sqrt(ui[0][i]*ui[0][i]+vi[0][i]*vi[0][i]
                             +wi[0][i]*wi[0][i]) / np.sqrt(params['gamma']*params['R'] * Tp))  

    # Plot Reynolds   
    #pl.plot_Mach(yi, Mach1p, params['path'])  

    # Print data
    for i in range(params['nfiles']):
        print('File'+str(i)+':')
        print('P2/P1=',P2[i]/P1[i])
        print('exit angle=', alphaM[i]*180/np.pi)
        print('mass avg exit angle=', alphaA[i]*180/np.pi)
        print('inlet angle=', alphaMI[i]*180/np.pi)
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
    #cp = fc.pressure_coefficient(p, P1, P2, npoints, params)

    # Axial chord
    x = objField.GetX()
    cax = fc.axial_chord(x[0]) # first file only

    # Get experimental data
    xexp = objExp.GetX()
    cpexp = objExp.GetCP()
    # Plot Cp distribution
    #pl.plot_cp(x/cax, cp, xexp, cpexp, params['path0']) 


    objExpCf = DataExpCf(params) # Experimental data
    # Get experimental data
    xexp = objExpCf.GetX()
    cfexp = objExpCf.GetCF()
    
    # X-shear stress 
    wss = objField.GetWSS()
    # Skin friction: Implementation from Garai et al. 2015
    cf = fc.skin_friction_coefficiet(x, wss, P1, P2, params) 
    # Plot x-shear stress
    #pl.plot_cf(x/cax, cf, xexp, cfexp, params['path0'])

    # rho density
    rho = objField.GetRho()
    # Calculate wall unit
    (spc, wallx, wally, wallz) = fc.wall_units_full(x, wss, rho, params) 
    # Plot wall units
    pl.plot_wall_unit(spc/cax, wallx, wally, wallz, params)

    # Wake loss
    ye = objExit.GetY() 
    npoints = objExit.GetNpoints()
    loss = fc.wake_loss(pe, rhoe, ue, ve, we, uM, P1, P2, npoints, params)
    ystar = fc.norm_position(ye, params['yminLoss'], params['Py'], params)

    # Exp
    objExpLoss = DataExpLoss(params) # Experimental data
    # Get experimental data
    xexp = objExpLoss.GetX()
    lossexp = objExpLoss.GetLoss()

    (lossShift, ystarShift) = fc.shift_data(loss, ystar, params['shiftLoss'])
    ystarFlip = fc.flip_data(ystarShift, params)

    # Plot wake loss
    #pl.plot_wake_loss(ystarFlip, lossShift, xexp, lossexp, params['path0'])

    # Mixed-out loss
    #mo_loss = fc.mixed_out_average_quantity(ystarFlip[0], -lossShift[0], params['Py'])
    #print("mo_loss = ", mo_loss)