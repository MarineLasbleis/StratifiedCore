#!/usr/bin/env python
## Time-stamp: "2015-12-03 17:03:20 marine"

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


IRON = {'name' : 'Fe', 'f1' : 5.26e-9, 'f2' : 1.24, 'f3' : -3.21}

def rho(f,f1,f2,f3) :
    """ equation 9 Gomi 2013 """
    return f1*(f2 - f)**f3
    


def resistivity(f, T, El, theta0=417., gamma=1.52, Ir=IRON) :
    """ Ideal resistivity from the Bloch GRuneisen formula as in Gomi 2013 (equation (1))
    
    f: V/V0
    T : absolute temperature
    El: element (dictionnary with z valence charge and x atomic fraction)
    theta0
    gamma
    
    """

    def integrand_eq1(z):
        """ integrande for function resistivity """
        return z**5/(np.exp(z)-1)/(1-np.exp(-z))

    theta = theta0*np.exp(-gamma *np.log(f))  #Debye temperature
    tfac = (T/theta)**5 *quad(integrand_eq1, 0, theta/T)[0]
    
    imp = El['z']**2*El['x']

    idSi = rho(f, El['f1'], El['f2'], El['f3'])
    idFe = rho(f, Ir['f1'], Ir['f2'], Ir['f3']) * tfac

    return idFe + idSi*imp  #ideal resistivity


def k(id, T, f, L=2.44e-8):
    """ real resistivity, adding ideal resistivity (id) and saturation resistivity (sat)
    saturation resistivity at 1bar is 1.68e-6 \Ohm m
    """
    sat = 1.68e-6 * f**(1./3.)
    return L*T*(1./id+ 1./sat)



if __name__ == '__main__':

    import eos

    ## Light elements :
    # name, z: valence charge, x: atomic fraction
    # f1,f2,f3 are fitting parameters
    IRON = {'name' : 'Fe', 'f1' : 5.26e-9, 'f2' : 1.24, 'f3' : -3.21}
    SILICIUM = {'name' : 'Si', 'z' : 1, 'x' : 22.5, 'f1' : 3.77e-8, 'f2' : 1.48, 'f3' : -3.10}
    CARBON = {'name' : 'C', 'z' : 1, 'x' : 30}
    OXYGEN = {'name' : 'O', 'z' : 0.5, 'x' : 23.2}
    SULFUR = {'name' : 'S', 'z' : 0.5, 'x' : 19.4}
    
    MURNAGHAN = {'T0': 1812., 'rho0': 7010., 'K': 130.*1e9, 'KPrim':4.,
                        'alphaPrim': 130.*np.log(2.)/(360.-135.) ,
                        'alpha': 1.e-5/np.exp(-0.400485 *135./130.),   #1.e-5*np.exp(130.*np.log(2.)/(360.-135.)*135./130),#   1.e-5/np.exp( np.log(2)/(360.-135.)*135.),
                        'Cp': 800., 'cond': 140.}
    P0 = 135e9
    T0 = 3750.

    N = 100 #number of points in the profile
    P = np.linspace(P0,330.e9,N)
    T = []

    ## Adiabatic profile
    for p in P:
        #print Tadiabatic(p, T0, P0, MURNAGHAN)
        T = np.append(T, eos.Tadiabatic(p, T0, P0, MURNAGHAN))

        
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(P/1.e9,T)
    ax1.set_xlabel("Pressure (GPa)")
    ax1.set_xlabel("Temperature (K)")
    ax1.set_title("Adiabatic profile")

    ## Associated k profile
    conductivity = []
    for i in range(0,N):        
        f = eos.fval(P[i],T[i], MURNAGHAN)
        res = resistivity(f, T[i], El=SILICIUM, theta0=417, gamma=1.52, Ir=IRON)
        conductivity = np.append(conductivity, k(res, T[i], f, L=2.44e-8))

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(P/1.e9,conductivity)
    ax2.set_xlabel("Pressure (GPa)")
    ax2.set_xlabel("Conductivity")
    ax2.set_title("Conductivity profile along the adiabatic profile")

    
    plt.show()
