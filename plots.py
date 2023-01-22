import matplotlib.pyplot as plt
from matplotlib import rc, cm
rc('text', usetex=True)

def plot_cp(x, cp, xe, cpe):
    # Figure - Cp distribution
    fig = plt.figure()
    fig.add_subplot(1,1,1)
    plt.rc('legend',**{'fontsize':14})

    p2, = plt.plot(xe[0], cpe[0], 'o', color='w', 
    markeredgecolor='k', markersize=11)	
    p0, = plt.plot(x[0], cp[0], '-', color='r', linewidth=3)	
    p1, = plt.plot(x[1], cp[1], '--', color='b', linewidth=3)	
    
    plt.legend([p0, p1, p2],
    [r'Nektar++: Stagnation Inflow bc',
    r'Nektar++: Enforce SU bc', 
    r'Experiment'],
    loc='best')
    
    plt.tick_params(reset=True, direction="in", which='both')
    plt.subplots_adjust(left=0.14, right=0.95, bottom=0.155, top=0.97)
    plt.grid(color='0.5', linestyle=':', linewidth=0.5, which='both')
    plt.xlim((0.0,1))
    plt.ylim((-0.5, 1.05))
    plt.xticks(fontsize = 24)
    plt.yticks(fontsize = 24)
    plt.xlabel(r'$x/C_{ax}$',fontsize = 24)
    plt.ylabel(r'$C_{p}$',fontsize = 24)
    fig.savefig('cp.pdf', format='PDF')
    fig.savefig('cp.png', format='png')
    plt.show()
