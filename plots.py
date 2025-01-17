import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, cm
rc('text', usetex=True)

def plot_tu_BL(s,tu,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':14})

    p0, = plt.plot(tu[0][0], s[0][0], '-', color='darkorange', linewidth=3)	
    plt.plot(tu[0][1], s[0][1], '-', color='darkorange', linewidth=3)	
    plt.plot(tu[0][2], s[0][2], '-', color='darkorange', linewidth=3)	
    plt.plot(tu[0][3], s[0][3], '-', color='darkorange', linewidth=3)	
    plt.plot(tu[0][4], s[0][4], '-', color='darkorange', linewidth=3)	
    plt.plot(tu[0][5], s[0][5], '-', color='darkorange', linewidth=3)	
    plt.plot(tu[0][6], s[0][6], '-', color='darkorange', linewidth=3)	
    plt.plot(tu[0][7], s[0][7], '-', color='darkorange', linewidth=3)

    p1, = plt.plot(tu[1][0], s[1][0], '-', color='m', linewidth=3)	
    plt.plot(tu[1][1], s[1][1], '-', color='m', linewidth=3)	
    plt.plot(tu[1][2], s[1][2], '-', color='m', linewidth=3)	
    plt.plot(tu[1][3], s[1][3], '-', color='m', linewidth=3)	
    plt.plot(tu[1][4], s[1][4], '-', color='m', linewidth=3)	
    plt.plot(tu[1][5], s[1][5], '-', color='m', linewidth=3)	
    plt.plot(tu[1][6], s[1][6], '-', color='m', linewidth=3)	
    plt.plot(tu[1][7], s[1][7], '-', color='m', linewidth=3)


    p2, = plt.plot(tu[2][0], s[2][0], '-', color='darkred', linewidth=3)	
    plt.plot(tu[2][1], s[2][1], '-', color='darkred', linewidth=3)	
    plt.plot(tu[2][2], s[2][2], '-', color='darkred', linewidth=3)	
    plt.plot(tu[2][3], s[2][3], '-', color='darkred', linewidth=3)	
    plt.plot(tu[2][4], s[2][4], '-', color='darkred', linewidth=3)	
    plt.plot(tu[2][5], s[2][5], '-', color='darkred', linewidth=3)	
    plt.plot(tu[2][6], s[2][6], '-', color='darkred', linewidth=3)
    plt.plot(tu[2][7], s[2][7], '-', color='darkred', linewidth=3)

    p3, = plt.plot(tu[3][0], s[3][0], '-', color='b', linewidth=3)	
    plt.plot(tu[3][1], s[3][1], '-', color='b', linewidth=3)	
    plt.plot(tu[3][2], s[3][2], '-', color='b', linewidth=3)	
    plt.plot(tu[3][3], s[3][3], '-', color='b', linewidth=3)	
    plt.plot(tu[3][4], s[3][4], '-', color='b', linewidth=3)	
    plt.plot(tu[3][5], s[3][5], '-', color='b', linewidth=3)	
    plt.plot(tu[3][6], s[3][6], '-', color='b', linewidth=3)
    plt.plot(tu[3][7], s[3][7], '-', color='b', linewidth=3)


    #plt.legend([p0,p1,p2,p3],
    #[
    #  r'0\% Tu',
    #  r'3.2\% Tu, 2\% chord length',
    #  r'3.2\% Tu, 5\% chord length',
    #  r'3.2\% Tu, 8\% chord length',
    #],
    #loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.152, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((0.55,1.175)) 
    plt.ylim((0.0, 0.0725)) 
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlabel(r'$x/C_{ax}$',fontsize = 24)
    plt.ylabel(r'$s/C$',fontsize = 24)
    fig.savefig(path + '/BLtu.pdf', format='PDF')
    fig.savefig(path + '/BLtu.png', format='png')
    plt.show()

def plot_up_BL(s,up,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':14})

    p0, = plt.plot(up[0][0], s[0][0], '-', color='darkorange', linewidth=3)	
    plt.plot(up[0][1], s[0][1], '-', color='darkorange', linewidth=3)	
    plt.plot(up[0][2], s[0][2], '-', color='darkorange', linewidth=3)	
    plt.plot(up[0][3], s[0][3], '-', color='darkorange', linewidth=3)	
    plt.plot(up[0][4], s[0][4], '-', color='darkorange', linewidth=3)	
    plt.plot(up[0][5], s[0][5], '-', color='darkorange', linewidth=3)	
    plt.plot(up[0][6], s[0][6], '-', color='darkorange', linewidth=3)	
    plt.plot(up[0][7], s[0][7], '-', color='darkorange', linewidth=3)

    p1, = plt.plot(up[1][0], s[1][0], '-', color='m', linewidth=3)	
    plt.plot(up[1][1], s[1][1], '-', color='m', linewidth=3)	
    plt.plot(up[1][2], s[1][2], '-', color='m', linewidth=3)	
    plt.plot(up[1][3], s[1][3], '-', color='m', linewidth=3)	
    plt.plot(up[1][4], s[1][4], '-', color='m', linewidth=3)	
    plt.plot(up[1][5], s[1][5], '-', color='m', linewidth=3)	
    plt.plot(up[1][6], s[1][6], '-', color='m', linewidth=3)	
    plt.plot(up[1][7], s[1][7], '-', color='m', linewidth=3)

    p2, = plt.plot(up[2][0], s[2][0], '-', color='darkred', linewidth=3)	
    plt.plot(up[2][1], s[2][1], '-', color='darkred', linewidth=3)	
    plt.plot(up[2][2], s[2][2], '-', color='darkred', linewidth=3)	
    plt.plot(up[2][3], s[2][3], '-', color='darkred', linewidth=3)	
    plt.plot(up[2][4], s[2][4], '-', color='darkred', linewidth=3)	
    plt.plot(up[2][5], s[2][5], '-', color='darkred', linewidth=3)	
    plt.plot(up[2][6], s[2][6], '-', color='darkred', linewidth=3)
    plt.plot(up[2][7], s[2][7], '-', color='darkred', linewidth=3)

    p3, = plt.plot(up[3][0], s[3][0], '-', color='b', linewidth=3)	
    plt.plot(up[3][1], s[3][1], '-', color='b', linewidth=3)	
    plt.plot(up[3][2], s[3][2], '-', color='b', linewidth=3)	
    plt.plot(up[3][3], s[3][3], '-', color='b', linewidth=3)	
    plt.plot(up[3][4], s[3][4], '-', color='b', linewidth=3)	
    plt.plot(up[3][5], s[3][5], '-', color='b', linewidth=3)	
    plt.plot(up[3][6], s[3][6], '-', color='b', linewidth=3)
    plt.plot(up[2][7], s[2][7], '-', color='b', linewidth=3)

    plt.legend([p0,p1,p2,p3],
    [
      r'0\% Tu',
      r'3.2\% Tu, 2\% chord length',
      r'3.2\% Tu, 5\% chord length',
      r'3.2\% Tu, 8\% chord length',
    ],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.152, right=0.96, bottom=0.135, top=0.96)
    plt.xlim((0.55,1.175)) 
    plt.ylim((0.0, 0.0725)) 
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlabel(r'$x/C_{ax}$',fontsize = 24)
    plt.ylabel(r'$s/C$',fontsize = 24)
    fig.savefig(path + '/BLup.pdf', format='PDF')
    fig.savefig(path + '/BLup.png', format='png')
    plt.show()

def plot_Tu(y,Tu,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':14})

    #p0, = plt.plot(Tu[0], y[0], '-', color='orangered', linewidth=3)	
    #p1, = plt.plot(Tu[1], y[0], '-', color='red', linewidth=3)	
    #p2, = plt.plot(Tu[2], y[0], '-', color='darkred', linewidth=3)	
    p3, = plt.plot(Tu[0], y[0], '-', color='darkred', linewidth=3)	


    plt.legend([p3],
    [
      #r'Present, 1.2\% Tu',
      #r'Present, 2\% Tu',
      #r'Present, 3.2\% Tu',
      r'Present, 3.8\% Tu',
    ],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.152, right=0.96, bottom=0.135, top=0.96)
    #plt.xlim((1e0,500)) 
    #plt.ylim((-150, 10)) 
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'Tu [$\%$]',fontsize = 18)
    plt.ylabel(r'$y$',fontsize = 18)
    fig.savefig(path + '/Tu.pdf', format='PDF')
    fig.savefig(path + '/Tu.png', format='png')
    plt.show()

def plot_tau(y,tau,yexp,tauexp,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':14})

    a = [0.0]
    b = [0.0]

    # Nektar++
    p0, = plt.plot(tau[0], y[0], '-', color='darkorange', linewidth=3)	
    p1, = plt.plot(tau[1], y[0], '-', color='red', linewidth=3)	
    p2, = plt.plot(tau[2], y[0], '-', color='crimson', linewidth=3)	
    p3, = plt.plot(tau[3], y[0], '-', color='darkred', linewidth=3)	
    # Exp
    plt.plot(tauexp[0], yexp[0], ':', color='darkorange', linewidth=3)	
    plt.plot(tauexp[1], yexp[1], ':', color='red', linewidth=3)	
    plt.plot(tauexp[2], yexp[2], ':', color='crimson', linewidth=3)	
    plt.plot(tauexp[3], yexp[3], ':', color='darkred', linewidth=3)	
    # Point
    p4, = plt.plot(a, b, ':', color='k', linewidth=3)	
     
    plt.legend([p4,p0,p1,p2,p3],
    [
     r'Sandberg $\&$ Michelassi 2019, 3.8\% Tu',
     r'$\tau_{11}$, 3.8\% Tu',
     r'$\tau_{22}$, 3.8\% Tu',
     r'$\tau_{33}$, 3.8\% Tu',
     r'$\tau_{12}$, 3.8\% Tu',
     ],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.152, right=0.96, bottom=0.135, top=0.96)
    #plt.xlim((1e0,500)) 
    #plt.ylim((-150, 10)) 
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$\tau_{ij}$',fontsize = 18)
    plt.ylabel(r'$y^{*}/P_{y}$',fontsize = 18)
    fig.savefig(path + '/tau.pdf', format='PDF')
    fig.savefig(path + '/tau.png', format='png')
    plt.show()

def plot_tke(y,tke,yexp,tkeexp,path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':14})

    #p1, = plt.plot(tke[0], y[0], '-', color='orangered', linewidth=3)	
    #p2, = plt.plot(tke[1], y[0], '-', color='red', linewidth=3)	
    #p3, = plt.plot(tke[2], y[0], '-', color='darkred', linewidth=3)	
    p4, = plt.plot(tke[0], y[0], '-', color='darkred', linewidth=3)	
    p0, = plt.plot(tkeexp[0], yexp[0], ':', color='k', linewidth=3)	

    plt.legend([p0,p4],
    [
      r'Sandberg $\&$ Michelassi 2019, 3.8\% Tu',
      #r'Present, 1.2\% Tu',
      #r'Present, 2\% Tu',
      #r'Present, 3.2\% Tu',
      r'Present, 3.8\% Tu',
    ],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.152, right=0.96, bottom=0.135, top=0.96)
    #plt.xlim((1e0,500)) 
    #plt.ylim((-150, 10)) 
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$TKE$',fontsize = 18)
    plt.ylabel(r'$y^{*}/P_{y}$',fontsize = 18)
    fig.savefig(path + '/tke.pdf', format='PDF')
    fig.savefig(path + '/tke.png', format='png')
    plt.show()

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
    fig.savefig(params['path'+str(params['nfiles']-1)] + '/psdu.pdf', format='PDF')
    fig.savefig(params['path'+str(params['nfiles']-1)] + '/psdu.png', format='png')
    plt.show()

def plot_psd_velocity_upstream(f,S,params):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':14})

    npsd = params['nfilepsdUp']

    p0, = plt.semilogx(f[npsd-4], S[npsd-4], '-', color='darkorange', linewidth=3)	
    p1, = plt.semilogx(f[npsd-3], S[npsd-3], '-', color='orangered', linewidth=3)	
    p2, = plt.semilogx(f[npsd-2], S[npsd-2], '-', color='red', linewidth=3)	
    p3, = plt.semilogx(f[npsd-1], S[npsd-1], '-', color='darkred', linewidth=3)	

    # -5/3 power law
    x0=0.3e1
    x1=2e2
    x = np.linspace(x0, x1)
    y = -(5/3)*10*np.log10(x) - (5/3)*10*np.log10(200) 
    p4, = plt.semilogx(x, y, '--', color='k')

    plt.legend([p0,p1,p2,p3,p4],
    [r'$-0.325C$',
     r'$-0.3C$',
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
    fig.savefig(params['path'+str(params['nfiles']-1)] + '/psdu.pdf', format='PDF')
    fig.savefig(params['path'+str(params['nfiles']-1)] + '/psdu.png', format='png')
    plt.show()

def plot_tpcorr(z, Rxx, Ryy, Rzz, path):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':20})

    #x = [-10]
    #y = [-10]
    # Fake plot
    #p0, = plt.plot(x, y, 'o-', color='k', markersize='6')
    #p1, = plt.plot(x, y, 's--', color='k', markersize='6')
    #p2, = plt.plot(x, y, '^--', color='k', markersize='6')
    #p3, = plt.plot(x, y, '-', color='m', markersize='6')
    #p4, = plt.plot(x, y, '-', color='darkred', markersize='6')
    #p5, = plt.plot(x, y, '-', color='b', markersize='6')

    # Real plot
    #p6, = plt.plot(z, Rxx[0], 'o-', color='m', markersize='6')	
    #p7, = plt.plot(z, Rxx[1], 'o-', color='darkred', markersize='6')	
    #p8, = plt.plot(z, Rxx[2], 'o-', color='b', markersize='6')	

    #p9, = plt.plot(z, Ryy[0], 's--', color='m', markersize='6')	
    #p10, = plt.plot(z, Ryy[1], 's--', color='darkred', markersize='6')	
    #p11, = plt.plot(z, Ryy[2], 's--', color='b', markersize='6')	

    #p12, = plt.plot(z, Rzz[0], '^:', color='m', markersize='6')	
    #p13, = plt.plot(z, Rzz[1], '^:', color='darkred', markersize='6')	
    #p14, = plt.plot(z, Rzz[2], '^:', color='b', markersize='6')	

    #plt.legend([p0, p1, p2],
    #[r'$R_{u^{\prime}u^{\prime}}$',
    # r'$R_{v^{\prime}v^{\prime}}$',
    # r'$R_{w^{\prime}w^{\prime}}$'],
    ## r'DNS data'],
    #loc='best')

    p0, = plt.plot(z, Rzz[0], '-', color='orange', linewidth='4')	
    p1, = plt.plot(z, Rzz[1], '-', color='m', linewidth='4')	
    p2, = plt.plot(z, Rzz[2], '-', color='darkred', linewidth='4')	
    p3, = plt.plot(z, Rzz[3], '-', color='b', linewidth='4')	

    #plt.legend([p0, p1, p2, p3],
    #[r'$0\%$ Tu',
    # r'$3.2\%$ Tu, $2\%$ chord length',
    # r'$3.2\%$ Tu, $5\%$ chord length',
    # r'$3.2\%$ Tu, $8\%$ chord length'],
    ## r'DNS data'],
    #loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='major')
    plt.subplots_adjust(left=0.148, right=0.98, bottom=0.17, top=0.98)
    #plt.xlim((-0.05,1.05))
    plt.ylim((-0.1, 1.05))
    plt.xticks(fontsize = 28)
    plt.yticks(fontsize = 28)
    plt.xlabel(r'$\Delta z / C$',fontsize = 28)
    plt.ylabel(r'$R_{w^{\prime}w^{\prime}}$',fontsize = 28)
    fig.savefig(path + '/corr.pdf', format='PDF')
    fig.savefig(path + '/corr.png', format='png')
    plt.show()

def plot_cp(x, cp, xe, cpe, path):
    # Figure - Cp distribution
    plt.figure()
    plt.rc('legend',**{'fontsize':12})

    #p0, = plt.plot(xe[0], cpe[0], 'o', color='none', markeredgecolor='k', markersize=13)	
    ##p1, = plt.plot(xe[1], cpe[1], ':', color='gold', linewidth=2)	
    #p2, = plt.plot(xe[2], cpe[2], '--', color='b', linewidth=2)	
    #p3, = plt.plot(xe[3], cpe[3], '--', color='purple', linewidth=2)	
    #p4, = plt.plot(x[0], cp[0], '-', color='darkorange', linewidth=2)	
    #p5, = plt.plot(x[1], cp[1], '-', color='orangered', linewidth=2)	
    #p6, = plt.plot(x[2], cp[2], '-', color='red', linewidth=2)	
    #p7, = plt.plot(x[3], cp[3], '-', color='darkred', linewidth=2)	

    p0, = plt.plot(xe[0], cpe[0], 'o', color='none', markeredgecolor='k', markersize=13)	
    ##p1, = plt.plot(xe[1], cpe[1], ':', color='gold', linewidth=2)	
    #p2, = plt.plot(xe[2], cpe[2], '--', color='b', linewidth=2)	
    #p3, = plt.plot(xe[3], cpe[3], '--', color='purple', linewidth=2)	
    p4, = plt.plot(x[0], cp[0], '-', color='darkorange', linewidth=2)	
    p5, = plt.plot(x[1], cp[1], '-', color='orangered', linewidth=2)	
    p6, = plt.plot(x[2], cp[2], '-', color='red', linewidth=2)	
    p7, = plt.plot(x[3], cp[3], '-', color='darkred', linewidth=2)	

    plt.legend([p0, p4, p5, p6, p7],
    [
     r'Experiment, 0\% Tu', 
    # #r'Wissink et al. 2003, 0\% Tu',
    # r'Michelassi et al. 2014, 0.0\% Tu',
    # r'Michelassi et al. 2014, 3.2\% Tu',
     r'Present, 0\% Tu',
     r'Present, 1.2\% Tu',
     r'Present, 2\% Tu',
     r'Present, 3.2\% Tu',
    ],
    loc='center')

    #p0, = plt.plot(xe[0], cpe[0], 'o', color='none', markeredgecolor='k', markersize=13)	
    #p1, = plt.plot(x[0], cp[0], '-', color='g', linewidth=2)	
    #p2, = plt.plot(x[1], cp[1], '-', color='darkred', linewidth=2)	
    #p3, = plt.plot(x[2], cp[2], '-', color='b', linewidth=2)	

    #plt.legend([p0, p1, p2, p3],
    #[
    # r'Experiment, 0\% Tu',
    # r'Present, 3.2\% Tu, 2\% chord length',
    # r'Present, 3.2\% Tu, 5\% chord length',
    # r'Present, 3.2\% Tu, 8\% chord length',
    #],
    #loc='center')
        
    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.18, right=0.95, bottom=0.16, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='both')
    plt.xlim((0.0,1))
    plt.ylim((-0.5, 1.05))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$x/C_{ax}$',fontsize = 20)
    plt.ylabel(r'$C_{p}$',fontsize = 20)
    plt.savefig(path + '/cp.pdf', format='PDF')
    plt.savefig(path + '/cp.png', format='png')
    plt.show()

def plot_cf(x, sh, xe, she, path):
    # Figure - Cp distribution
    plt.figure()
    plt.rc('legend',**{'fontsize':16})

    # x-axis: x=0 line
    a = [0, 1]
    b = [0, 0]
    plt.plot(a, b, ':', color='k')

    # Just suction side
    shaux = []
    xaux = []
    for i in range(len(sh)):
        shaux.append([])
        xaux.append([])
        for j in range(400, len(sh[i])):
            shaux[i].append(sh[i][j])
            xaux[i].append(x[i][j])

    for i in range(len(sh)):
        for j in range(0, 147):
            shaux[i].append(sh[i][j])
            xaux[i].append(x[i][j])


    #p0, = plt.plot(xe[0], she[0], '--', color='purple', linewidth=2)
    #p0, = plt.plot(xaux[0], shaux[0], '-', color='darkorange', linewidth=2)	
    #p1, = plt.plot(xaux[1], shaux[1], '-', color='orangered', linewidth=2)	
    p0, = plt.plot(xaux[0], shaux[0], '-', color='m', linewidth=2)	
    p1, = plt.plot(xaux[1], shaux[1], '-', color='darkred', linewidth=2)	
    p2, = plt.plot(xaux[2], shaux[2], '-', color='b', linewidth=2)	

    plt.legend([p0, p1, p2],
    [
     r'3.2\% Tu, 2\% chord length',
     r'3.2\% Tu, 5\% chord length',
     r'3.2\% Tu, 8\% chord length',
    ],
    loc='best')


    #plt.legend([p0, p1, p2, p3],
    #[
    # r'0\% Tu',
    # r'1.2\% Tu',
    # r'2\% Tu',
    # r'3.2\% Tu',
    # ],
    #loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.245, right=0.95, bottom=0.16, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='both')
    plt.xlim((0.6,1))
    plt.ylim((-0.001, 0.004))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$x/C_{ax}$',fontsize = 20)
    plt.ylabel(r'$C_{f}$',fontsize = 20)
    plt.savefig(path + '/cf.pdf', format='PDF')
    plt.savefig(path + '/cf.png', format='png')
    plt.show()

def plot_wake_loss(x, loss, xe, losse, path):
    # Figure - Cp distribution
    plt.figure()
    plt.rc('legend',**{'fontsize':12})

    p0, = plt.plot(xe[0], losse[0], 'o', color='none', 
                   markeredgecolor='k', markersize=13)	
    #p1, = plt.plot(xe[1], losse[1], '--', color='b', linewidth=2)	#0% TI M
    ##p2, = plt.plot(xe[3], losse[3], '--', color='m', linewidth=2)	#1.2% TI M
    #p3, = plt.plot(xe[4], losse[4], '--', color='purple', linewidth=2)	#3.2% TI M

    ##p4, = plt.plot(xe[2], losse[2], ':', color='g', linewidth=2)  #0% TI G	
    ##p5, = plt.plot(xe[5], losse[5], ':', color='springgreen', linewidth=2)  #0% TI G	

    #p6, = plt.plot(x[0], loss[0], '-', color='darkorange', linewidth=2)	
    #p7, = plt.plot(x[1], loss[1], '-', color='orangered', linewidth=2)	
    #p8, = plt.plot(x[2], loss[2], '-', color='red', linewidth=2)	
    #p9, = plt.plot(x[3], loss[3], '-', color='darkred', linewidth=2)	

    # Different length scales
    p1, = plt.plot(x[0], loss[0], '-', color='m', linewidth=2)	
    p2, = plt.plot(x[1], loss[1], '-', color='darkred', linewidth=2)	
    p3, = plt.plot(x[2], loss[2], '-', color='b', linewidth=2)	


#    plt.legend([p0, p1, p3, p6, p7,p8, p9],
#    [
#     r'Experiment, 0\% Tu',
#     r'Michelassi et al. 2014, 0\% Tu',
#     #r'Michelassi et al. 2014, 1.2\% Tu',
#     r'Michelassi et al. 2014, 3.2\% Tu',
#     #r'Garai et al. 2016, 0\% Tu',
#     #r'Garai et al. 2016, 2\% Tu',
#     r'Present, 0\% Tu',
#     r'Present, 1.2\% Tu',
#     r'Present, 2\% Tu',
#     r'Present, 3.2\% Tu',
#    ],
#    loc='best')

    plt.legend([p0, p1, p2, p3],
    [
     r'Experiment, 0\% Tu',
     r'Present, 3.2\% Tu, 2\% chord length',
     r'Present, 3.2\% Tu, 5\% chord length',
     r'Present, 3.2\% Tu, 8\% chord length',
    ],
    loc='best')

    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.16, right=0.95, bottom=0.16, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='both')
    plt.xlim((-0.02,1.02))
    plt.ylim((-0.008, 0.25))
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$y^{*}/P_{y}$',fontsize = 20)
    plt.ylabel(r'$\Omega$',fontsize = 20)
    plt.savefig(path + '/loss.pdf', format='PDF')
    plt.savefig(path + '/loss.png', format='png')
    plt.show()

def plot_wall_unit(spc, wallx, wally, wallz, params):
    # Figure - Cp distribution
    plt.figure()
    plt.rc('legend',**{'fontsize':18})

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
    plt.text(0.06, 22, r'$\Delta x^{+}$ DNS limit', color='darkorange',fontsize = 16)
    plt.text(0.06, 1.13, r'$\Delta y^{+}$ DNS limit', color='red',fontsize = 16)
    plt.text(0.06, 11, r'$\Delta z^{+}$ DNS limit', color='darkred',fontsize = 16)


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
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.xlabel(r'$x/C_{ax}$',fontsize = 20)
    plt.ylabel(r'Wall resolution',fontsize = 20)
    plt.savefig(params['path0'] + '/wu.pdf', format='PDF')
    plt.savefig(params['path0'] + '/wu.png', format='png')
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