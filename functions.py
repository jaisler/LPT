import math
import numpy as np
from scipy import signal
import csv
import pandas as pd

def bubble_sort(y, qt1, qt2, qt3, qt4, qt5):
    "Bubble sort algorithm"
    has_swapped = True
    num_of_iterations = 0

    # Reverse array
    yr = y[::-1]
    qt1r = qt1[::-1]
    qt2r = qt2[::-1]
    qt3r = qt3[::-1]
    qt4r = qt4[::-1]
    qt5r = qt5[::-1]

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
                qt5r[i], qt5r[i+1] = qt5r[i+1], qt5r[i]
                has_swapped = True
        num_of_iterations += 1
    return(yr, qt1r, qt2r, qt3r, qt4r, qt5r)

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
    #for i in range(1, int(len(y)/2)): #Simpson's 1/3 Rule
    for i in range(len(y)-1):
        #num = num + ((qt[2*i-2]*u[2*i-2]+4*qt[2*i-1]*u[2*i-1]+qt[2*i]*u[2*i])*(y[2*i-2]-y[2*i]))/6.0
        #den = den + ((u[2*i-2]+4*u[2*i-1]+u[2*i])*(y[2*i-2]-y[2*i]))/6.0
        num += ((qt[i+1]*u[i+1]+qt[i]*u[i])*(y[i+1]-y[i]))/2.0
        den += ((u[i+1]+u[i])*(y[i+1]-y[i]))/2.0
    return(num/den)

def mixed_out_average_quantity(y, qt, pitch):
    """ Mixed-out pressure reference at 0.5*Cax downstream of the blade trailing edge """
    qtbar = 0.0
    for i in range(len(y)-1):
        qtbar += ((qt[i+1]+qt[i])*(y[i+1]-y[i]))/2.0
    qtbar = qtbar/pitch
    return(qtbar)

def pressure_coefficient(p, P1, P2, npoints, params):
    """ Calculate the pressure distribution on a blade """
    cp = []
    for j in range(params['nfiles']):
        cp.append([])
        for i in range(npoints[j]):
            cp[j].append((p[j][i]-(P2[j]))/(P1[j]-(P2[j]))) #P2+0.004
    return (cp)

def skin_friction_coefficiet(x, wss, P1, P2, params):
    """ Calculate the skin friction coefficient on a blade"""
    cf = []
    for i in range(params['nfiles']):
        cf.append([])
        flag1 = True
        flag2 = False
        sign = 1
        for j in range(len(wss[i]) - 1):
            if x[i][j] > x[i][j+1] and flag1:
                sign *= -1 
                flag1 = False
                flag2 = True
            elif x[i][j] < 0 and flag2:
                sign *= -1 
                flag2 = False
            cf[i].append(sign * wss[i][j]/(P1[i]-P2[i])) # sign
        cf[i].append((sign * wss[i][len(x[i])-1])/(P1[i]-P2[i]))

    return cf

def wake_loss(p, rhoe, ue, ve, we, uM, P1, P2, npoints, params):
    """ Calculate the wake loss for the blade """
    loss = []
    for j in range(params['nfiles']):
        loss.append([])
        for i in range(params['spline']):
            mod2 = 0.5 * rhoe[j][i] * (ue[j][i] * ue[j][i] + ve[j][i] * ve[j][i] + 
                                       we[j][i] * we[j][i])
            loss[j].append((P1[j] - 0.052 - (p[j][i] + mod2))/(0.5 * rhoe[j][i] * uM[j] * uM[j]))

    return (loss)

def shift_data(u, y, shiftPar):
    """ Calculate shift for data """
    shift = []
    for i in range(len(u)):
        shift.append(shiftPar)
    uShift = []
    yeShift = []
    for j in range(len(u)):
        uShift.append(np.zeros(len(u[j])))
        value = y[j][shift[j]] 
        k = shift[j]
        yeShift.append(np.zeros(len(u[j])))
        for i in range(len(u[j])):
            uShift[j][i] = u[j][k]
            yeShift[j][i] = y[j][k] - value
            if (k == (len(u[j]) - 1)):
                k = 0
                shift[j] = len(u[j]) - shift[j] 
                value = - y[j][shift[j]] #+ add[j]
            else:
                k += 1
    
    return (uShift, yeShift)

def flip_data(y, params):
    """ Flip data  """
    yFlip = []
    for j in range(params['nfiles']):
        k = len(y[j])
        yFlip.append([])
        for i in range(k):
            if i == 0:
                yFlip[j].append(y[j][k-i-1])
            else:
                yFlip[j].append(yFlip[j][i-1] - abs(y[j][i]-y[j][i-1]))

    return yFlip

def norm_position(y, ymin, Py, params):
    """ Calculate the correct position and then normalise the data """   
    ystar = []
    for j in range(len(y)):
        ystar.append((y[j] - ymin) / Py)
    return (ystar)


def CalculatePSD(u, dt, params, flag):
    """ Calculate PSD """
    # f contains the frequency components
    # S is the PSD
    pf = []
    pS = []
    fs = []

    if flag == 'down':
        nfiles = params['nfilepsd']
    elif flag == 'up':
        nfiles = params['nfilepsdUp']

    #print(math.ceil(params['nfiles']/2))
    # Obtain the time step for each file and adjust parameters
    for i in range(nfiles):
        dt[i] *= params['C']/params['Uinf']
        fs.append(int(1/dt[i]))
        u[i] *= params['Uinf']

    # Calculare PSD
    nameWin = 'hanning'
    for i in range(nfiles):
        nx = int(len(u[i])/3)  # /6 number of windows - 3 or 6
        (f, S) = signal.welch(u[i], fs[i], scaling='density', 
                                    window=signal.get_window(nameWin, nx),
                                    nperseg=nx, noverlap=round(nx/2))
        pf.append(f)
        pS.append(10.0*np.log10(S)) #  
        del f, S
    return (pf, pS)

def CalculateTwoPointsCorrelation(u, params):
    """ Calculate two-points correlation function """
    # denominator - RMS
    k = 0
    tn = 0
    den = 0
    while (k < len(u)):
        meanz = 0
        for i in range(params['npcorr']): 
            meanz += u[i + k] * u[i + k]
        meanz /= (params['npcorr'])
        den += meanz
        k += params['npcorr']
        tn += 1
    den /= tn

    # numerator
    Rii = np.zeros(params['npcorr']) # correlation
    for l in range(params['npcorr']):
        k = 0
        tn = 0
        while (k < len(u)): # time
            for i in range(params['npcorr']): # space
                if ((i + l) < params['npcorr']):
                    Rii[l] += u[i + k] * u[i + l + k] 
                else:
                    Rii[l] += u[i + k] * u[i + l + 2 * (params['npcorr'] - 1 - i - l) + k]
            k += params['npcorr']
            tn += 1
        Rii[l] /= tn
        Rii[l] /= (params['npcorr'] * den) 

    return Rii

def wall_units(x, wss, rho, params):
    """ Calculate wall units """
    wallx = []
    wally = []
    wallz = []
    spc = []
    shear = []
    den = []
    kinvis = []
    xi = 573 # by hand
    xf = 147 # by hand
    for i in range(params['nfiles']):
        spc.append([])
        shear.append([])
        kinvis.append([])
        den.append([])
        for j in range(xi, len(x[i])):
            spc[i].append(x[i][j])
            shear[i].append(wss[i][j])
            den[i].append(rho[i][j])
            kinvis[i].append(params['mu'][0]/rho[i][j])
        for j in range(xf):
            spc[i].append(x[i][j])
            shear[i].append(wss[i][j])
            den[i].append(rho[i][j])
            kinvis[i].append(params['mu'][0]/rho[i][j])
    
    # Calculation of kinvis with rhoInf
    kinvisStat = 1 # params['rhoInf'] / params['mu'][0]

    for i in range(params['nfiles']):
        wallx.append([])
        wally.append([])
        wallz.append([])
        for j in range(len(spc[0]) - 1):
            dx = spc[i][j+1] - spc[i][j]
            wallx[i].append(dx * np.sqrt(abs(shear[i][j])/rho[i][j]) / kinvis[i][j])
            wally[i].append((params['dyplus'] * np.sqrt(abs(shear[i][j])/rho[i][j])) / kinvis[i][j])
            wallz[i].append((params['dzplus'] * np.sqrt(abs(shear[i][j])/rho[i][j])) / kinvis[i][j])
        spc[i].pop()

    return (spc, wallx, wally, wallz)

def boundary_layer_adjustment(factor, u, params):
    """ Adjust boundary layer turbulence intensity """

    post = 0.50 # first point
    dec = 0.05
    for i in range(params['nfilesBLtu']):
        u[i] /= factor    
        u[i] += (post/params['Cax'])
        post += 0.05

    return (u)


def WriteDataToCSV(filename, headers, s, up, vp, ump, taumz, tke, tu):
    """
    Writes two vectors to a CSV file with headers.

    Args:
    filename (str): The name of the file to write to.
    headers (list of str): The headers for the CSV file.
    vector (list): The vector (column).
    matrix (list): The matrix (columns).
    """

    # Open a file in write mode
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)

        # Write the header
        writer.writerow(headers)

        # Write the data
        for v0, v1, v2, v3, v4, v5, v6, v7, v8, v9 in zip(s, up, vp, ump, taumz[0], taumz[1], taumz[2], taumz[3], tke, tu):
            writer.writerow([v0,v1,v2,v3,v4,v5,v6,v7,v8,v9])

def LoadCSV(filename, delimiter=',', header='infer', skiprows=None, usecols=None):
    """
    Load data from a CSV file into a DataFrame.

    Args:
    filename (str): The path to the CSV file to load.
    delimiter (str, optional): The character used to separate fields. Defaults to ','.
    header (int, list of int, or str, optional): Row(s) to use as column names, and the start of data.
        Defaults to 'infer' which means pandas tries to guess the header location.
    skiprows (list-like, int or callable, optional): Line numbers to skip (0-indexed) or number of lines to skip (int) at the start of the file.
    usecols (list-like or callable, optional): Return a subset of the columns by specifying column names or indices.

    Returns:
    pandas.DataFrame: A DataFrame containing the loaded data.
    """
    return pd.read_csv(
        filename,
        delimiter=delimiter,
        header=header,
        skiprows=skiprows,
        usecols=usecols
    )