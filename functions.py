import math

def bubble_sort(y, qt1, qt2, qt3, qt4):
    "Bubble sort algorithm"
    has_swapped = True
    num_of_iterations = 0

    # Reverse array
    yr = y[::-1]
    qt1r = qt1[::-1]
    qt2r = qt2[::-1]
    qt3r = qt3[::-1]
    qt4r = qt4[::-1]

    while(has_swapped):
        has_swapped = False
        for i in range(len(yr) - num_of_iterations - 1):
            if yr[i] > yr[i+1]:
                # Swap
                yr[i], yr[i+1] = yr[i+1], yr[i]
                qt1r[i], qt1r[i+1] = qt1r[i+1], qt1r[i]
                qt2r[i], qt2r[i+1] = qt2r[i+1], qt2r[i]
                qt3r[i], qt3r[i+1] = qt3r[i+1], qt3r[i]
                qt4r[i], qt4r[i+1] = qt4r[i+1], qt4r[i]
                has_swapped = True
        num_of_iterations += 1
    return(yr, qt1r, qt2r, qt3r, qt4r)

def axial_chord(x):
    """ Obtain the axial chord of the blade """
    CaxMax = 0
    CaxMin = 0
    for i in range(len(x)):
        if CaxMax < x[i]:
            CaxMax = x[i]
        if CaxMin > x[i]:
            CaxMin = x[i]
    return(CaxMax-CaxMin) 

def arithmetic_mean(qt):
    "Arithmetic_mean"
    soma = 0
    for i in range(len(qt)):
        soma = soma + qt[i]
    return(soma/len(qt))

def mass_flow_average_quantity(y, u, qt):
    """ Mass average quantity """
    num = 0
    den = 0
    for i in range(1, int(len(y)/2)): #Simpson's 1/3 Rule
        num = num + ((qt[2*i-2]*u[2*i-2]+4*qt[2*i-1]*u[2*i-1]+qt[2*i]*u[2*i])*(y[2*i-2]-y[2*i]))/6.0
        den = den + ((u[2*i-2]+4*u[2*i-1]+u[2*i])*(y[2*i-2]-y[2*i]))/6.0
    return(num/den)

def mixed_out_average_quantity(y, qt, pitch):
    """ Mixed-out pressure reference at 0.5*Cax downstream of the blade trailing edge """
    qtbar = 0.0
    for i in range(len(y)-1):
        qtbar += ((qt[i+1]+qt[i])*(y[i+1]-y[i]))/2.0
    qtbar = qtbar/pitch
    return(qtbar)

def cp_distribution(p, P1, P2, npoints, params):
    """ Calculate the pressure distribution on a blade """
    cp = []
    for j in range(params['nfiles']):
        cp.append([])
        for i in range(npoints[j]):
            cp[j].append((p[j][i]-P2)/(P1-P2))
    return (cp)





def exit_Mach_number(y, u, v, pitch, cInf):
    "Obtain the exit Mach number at 0.5*Cax downstram of the blade trailing edge."
    (ys, us, vs) = bubble_sort(y, u, v)
    uM = mixed_out_velocity(ys, us, vs, pitch)
    return(uM/cInf)

def exit_Reynolds_number(y, u, v, rho, mu, pitch, C):
    #"Obtain the exit Mach number at 0.5*Cax downstram of the blade trailing edge."
    (ys, us, vs) = bubble_sort(y, u, v)
    (ys, us, rhos) = bubble_sort(y, u, rho) #Repetir por cause do rho
    uM = mixed_out_velocity(ys, us, vs, pitch)
    rhobar = mixed_out_density(ys, rhos, pitch)
    return((rhobar*C*uM)/mu)

def mixed_out_velocity(y, u, v, pitch):
    "Mixed-out velocity reference at 0.5*Cax downstream of the blade trailing edge."
    ubar = 0
    vbar = 0
    for i in range(len(y)-1):
        ubar = ubar + ((u[i+1]+u[i])*(y[i+1]-y[i]))/2.0
        vbar = vbar + ((v[i+1]+v[i])*(y[i+1]-y[i]))/2.0
    ubar = ubar/pitch
    vbar = vbar/pitch
    return(math.sqrt(ubar*ubar+vbar*vbar))

def mixed_out_density(y, rho, pitch):
    "Mixed-out density reference at 0.5*Cax downstream of the blade trailing edge."
    rhobar = 0
    for i in range(len(y)-1):
        rhobar = rhobar + ((rho[i+1]+rho[i])*(y[i+1]-y[i]))/2.0
    rhobar = rhobar/pitch
    return(rhobar)


