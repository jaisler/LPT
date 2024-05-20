import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, cm
rc('text', usetex=True)

def plot_psd_velocity(f,S,params):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':14})

    npsd = params['nfilepsd']

    #p0, = plt.semilogx(f[0], S[0], '-', color='darkorange', linewidth=3)	
    #p0, = plt.semilogx(f[1], S[1], '-', color='m', linewidth=3)	
    p0, = plt.semilogx(f[npsd-13], S[npsd-13], '-', color='darkorange', linewidth=3)	
    p1, = plt.semilogx(f[npsd-12], S[npsd-12], '-', color='orangered', linewidth=3)	
    p2, = plt.semilogx(f[npsd-3], S[npsd-3], '-', color='red', linewidth=3)	
    p3, = plt.semilogx(f[npsd-2], S[npsd-2], '-', color='crimson', linewidth=3)	
    p4, = plt.semilogx(f[npsd-1], S[npsd-1], '-', color='darkred', linewidth=3)	

    # -5/3 power law
    x0=0.3e1
    x1=2e2
    x = np.linspace(x0, x1)
    y = -(5/3)*10*np.log10(x) - (5/3)*10*np.log10(4) 
    p5, = plt.semilogx(x, y, '--', color='k')

    plt.legend([p0,p1,p2,p3,p4,p5],
    [r'Bubble $1$',
     r'Bubble $2$',
     r'Wake $1$',
     r'Wake $2$',
     r'Wake $3$',
     r'$f^{-5/3}$'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.152, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((1e0,500)) # 0.7,1000
    plt.ylim((-150, 10)) # -220, 20
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$f C / U_{\infty}$',fontsize = 18)
    plt.ylabel(r'$PSD(u)$ $[dB]$',fontsize = 18)
    fig.savefig(params['path'] + '/psdu.pdf', format='PDF')
    fig.savefig(params['path'] + '/psdu.png', format='png')
    plt.show()

def plot_psd_velocity_upstream(f,S,params):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':14})

    npsd = params['nfilepsdUp']

    p0, = plt.semilogx(f[npsd-3], S[npsd-3], '-', color='darkorange', linewidth=3)	
    p1, = plt.semilogx(f[npsd-2], S[npsd-2], '-', color='red', linewidth=3)	
    p2, = plt.semilogx(f[npsd-1], S[npsd-1], '-', color='darkred', linewidth=3)	

    # -5/3 power law
    x0=0.3e1
    x1=2e2
    x = np.linspace(x0, x1)
    y = -(5/3)*10*np.log10(x) - (5/3)*10*np.log10(200) 
    p3, = plt.semilogx(x, y, '--', color='k')

    plt.legend([p0,p1,p2,p3],
    [r'$-0.3C$',
     r'$-0.2C$',
     r'$-0.1C$',
     r'$f^{-5/3}$'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.152, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((1e0,500)) # 0.7,1000
    plt.ylim((-150, 10)) # -220, 20
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$f C / U_{\infty}$',fontsize = 18)
    plt.ylabel(r'$PSD(u)$ $[dB]$',fontsize = 18)
    fig.savefig(params['path'] + '/psdu.pdf', format='PDF')
    fig.savefig(params['path'] + '/psdu.png', format='png')
    plt.show()

def plot_tpcorr(z, R, path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':18})

    p0, = plt.plot(z, R[0], 'o-', color='r', markersize='6')	
    p1, = plt.plot(z, R[1], 's-', color='b', markersize='6')	
    p2, = plt.plot(z, R[2], '^-', color='m', markersize='6')	

    plt.legend([p0, p1, p2],
    [r'$R_{u^{\prime}u^{\prime}}$',
     r'$R_{v^{\prime}v^{\prime}}$',
     r'$R_{w^{\prime}w^{\prime}}$'],
    # r'DNS data'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.135, top=0.96)
    #plt.xlim((0,6.3))
    #plt.ylim((-0.1, 1.1))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$\Delta z / C$',fontsize = 18)
    plt.ylabel(r'$R_{i^{\prime}i^{\prime}}$',fontsize = 18)
    fig.savefig(path + '/corr.pdf', format='PDF')
    fig.savefig(path + '/corr.png', format='png')
    plt.show()

def plot_cp(x, cp, xe, cpe, path):
    # Figure - Cp distribution
    plt.figure()
    plt.rc('legend',**{'fontsize':12})

    p0, = plt.plot(xe[0], cpe[0], 'o', color='none', markeredgecolor='k', markersize=13)	
    p1, = plt.plot(xe[1], cpe[1], ':', color='gold', linewidth=2)	
    p2, = plt.plot(xe[2], cpe[2], '--', color='b', linewidth=2)	
    p3, = plt.plot(xe[3], cpe[3], '--', color='purple', linewidth=2)	
    p4, = plt.plot(x[0], cp[0], '-', color='red', linewidth=2)	
    p5, = plt.plot(x[1], cp[1], '-', color='darkorange', linewidth=2)	

    plt.legend([p0, p1, p2, p3, p4, p5],
    [r'Experiment, 0\% Tu', 
     r'Wissink et al. 2003, 0\% Tu',
     r'Michelassi et al. 2014, 0.0\% Tu',
     r'Michelassi et al. 2014, 3.2\% Tu',
     r'Present, 0\% Tu',
     r'Present, 3.2\% Tu'],
    loc='center')
    
    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.18, right=0.95, bottom=0.16, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='both')
    plt.xlim((0.0,1))
    plt.ylim((-0.5, 1.05))
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlabel(r'$x/C_{ax}$',fontsize = 24)
    plt.ylabel(r'$C_{p}$',fontsize = 24)
    plt.savefig(path + '/cp.pdf', format='PDF')
    plt.savefig(path + '/cp.png', format='png')
    plt.show()

def plot_cf(x, sh, xe, she, path):
    # Figure - Cp distribution
    plt.figure()
    plt.rc('legend',**{'fontsize':15})

    # x-axis: x=0 line
    a = [0, 1]
    b = [0, 0]
    plt.plot(a, b, ':', color='k')

    #p0, = plt.plot(xe[0], she[0], '--', color='purple', linewidth=2)	
    p0, = plt.plot(x[0], sh[0], '-', color='red', linewidth=2)	
    p1, = plt.plot(x[1], sh[1], '-', color='darkorange', linewidth=2)	

    plt.legend([p0, p1],
    #[r'Nektar++: Stagnation Inflow bc',
    #[r'Michelassi et al. 2014, 3.2\% Tu',
    [r'Present, 0\% Tu',
    r'Present, 3.2\% Tu'],
    #r'Present, 7.8\% Tu'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.245, right=0.95, bottom=0.16, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='both')
    plt.xlim((0.0,1))
    plt.ylim((-0.01, 0.01))
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlabel(r'$x/C_{ax}$',fontsize = 24)
    plt.ylabel(r'$C_{f}$',fontsize = 24)
    plt.savefig(path + '/cf.pdf', format='PDF')
    plt.savefig(path + '/cf.png', format='png')
    plt.show()

def plot_wake_loss(x, loss, xe, losse, path):
    # Figure - Cp distribution
    plt.figure()
    plt.rc('legend',**{'fontsize':10})

    p0, = plt.plot(xe[0], losse[0], 'o', color='none', 
                   markeredgecolor='k', markersize=13)	
    p1, = plt.plot(xe[1], losse[1], '--', color='b', linewidth=2)	#0% TI M
    p2, = plt.plot(xe[3], losse[3], '--', color='m', linewidth=2)	#1.2% TI M
    p3, = plt.plot(xe[4], losse[4], '--', color='purple', linewidth=2)	#3.2% TI M

    p4, = plt.plot(xe[2], losse[2], ':', color='g', linewidth=2)  #0% TI G	
    p5, = plt.plot(xe[5], losse[5], ':', color='springgreen', linewidth=2)  #0% TI G	

    p6, = plt.plot(x[0], loss[0], '-', color='red', linewidth=2)	
    p7, = plt.plot(x[1], loss[1], '-', color='darkorange', linewidth=2)	
    #p8, = plt.plot(x[2], loss[2], '-', color='darkorange', linewidth=2)	


    plt.legend([p0, p1, p2, p3, p4, p5, p6, p7],#, ],
    #[r'Nektar++: Stagnation Inflow bc',
    [r'Experiment, 0\% Tu',
    r'Michelassi et al. 2014, 0\% Tu',
    r'Michelassi et al. 2014, 1.2\% Tu',
    r'Michelassi et al. 2014, 3.2\% Tu',
    r'Garai et al. 2016, 0\% Tu',
    r'Garai et al. 2016, 2\% Tu',
    r'Present, 0\% Tu',
    r'Present, 3.2\% Tu'],
    #r'Present, 7.8\% Tu'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.16, right=0.95, bottom=0.16, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='both')
    plt.xlim((-0.02,1.02))
    plt.ylim((-0.008, 0.25))
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlabel(r'$y^{*}/P_{y}$',fontsize = 24)
    plt.ylabel(r'$\Omega$',fontsize = 24)
    plt.savefig(path + '/loss.pdf', format='PDF')
    plt.savefig(path + '/loss.png', format='png')
    plt.show()

def plot_wall_unit(spc, wallx, wally, wallz, params):
    # Figure - Cp distribution
    plt.figure()
    plt.rc('legend',**{'fontsize':15})

    p0, = plt.semilogy(spc[0], wallx[0], '--', color='darkorange', linewidth=2)	
    p1, = plt.semilogy(spc[0], wally[0], '-.', color='red', linewidth=2)	
    p2, = plt.semilogy(spc[0], wallz[0], ':', color='darkred', linewidth=2)	
    #p2, = plt.plot(xe[0], losse[0], 'o', color='none', 
    #               markeredgecolor='k', markersize=13)	

    x = [0,1]
    # Delta x plus
    xdns = [20,20]
    # Delta y plus
    ydns = [1,1]
    # Delta z plus
    zdns = [10,10]
    pa, = plt.semilogy(x, xdns, '--', color='darkorange', linewidth=1)	
    pb, = plt.semilogy(x, ydns, '-.', color='red', linewidth=1)	
    pc, = plt.semilogy(x, zdns, ':', color='darkred', linewidth=1)	

    # Test
    plt.text(0.06, 22, r'$\Delta x^{+}$ DNS limit', color='darkorange')
    plt.text(0.06, 1.13, r'$\Delta y^{+}$ DNS limit', color='red')
    plt.text(0.06, 11, r'$\Delta z^{+}$ DNS limit', color='darkred')


    plt.legend([p0, p1, p2],
    [r'$\Delta x^{+}$',
     r'$\Delta y_{wall}^{+}$', 
     r'$\Delta z^{+}$'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.16, right=0.95, bottom=0.16, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.xlim((0.0,1.0))
    plt.ylim((0.001, 50))
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlabel(r'$x/C_{ax}$',fontsize = 24)
    plt.ylabel(r'Wall units',fontsize = 24)
    plt.savefig(params['path'] + '/wu.pdf', format='PDF')
    plt.savefig(params['path'] + '/wu.png', format='png')
    plt.show()


def plot_Re(ye, Re, path):
    # Figure - exit isentropic Reynolds number
    plt.figure()
    plt.rc('legend',**{'fontsize':15})

    p0, = plt.plot(Re, ye[0], '-', color='r', linewidth=2)	
    #p1, = plt.plot(x[1], sh[1], '--', color='b', linewidth=2)	

    plt.legend([p0],
    [r'3D, 0\% TI'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.18, right=0.95, bottom=0.16, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='both')
    #plt.xlim((0.0,1))
    #plt.ylim((-0.01, 0.01))
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlabel(r'$Re_{2}^{is}$',fontsize = 24)
    plt.ylabel(r'$y$',fontsize = 24)
    plt.savefig(path + '/re.pdf', format='PDF')
    plt.savefig(path + '/re.png', format='png')
    plt.show()

def plot_Mach(ye, Me, path):
    # Figure - exit isentropic Reynolds number
    plt.figure()
    plt.rc('legend',**{'fontsize':15})

    p0, = plt.plot(Me, ye[0], '-', color='r', linewidth=2)	
    #p1, = plt.plot(x[1], sh[1], '--', color='b', linewidth=2)	

    plt.legend([p0],
    #[r'Nektar++: Stagnation Inflow bc',
    [r'3D, 0\% TI'],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.18, right=0.95, bottom=0.16, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='both')
    #plt.xlim((0.0,1))
    #plt.ylim((-0.01, 0.01))
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlabel(r'$Ma_{2}^{is}$',fontsize = 24)
    plt.ylabel(r'$y$',fontsize = 24)
    plt.savefig(path + '/mach.pdf', format='PDF')
    plt.savefig(path + '/mach.png', format='png')
    plt.show()