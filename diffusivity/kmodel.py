#!/usr/bin/env python

import numpy as np
from scipy.integrate import quad

## Light elements :
# name, z: valence charge, x: atomic fraction
# f1,f2,f3 are fitting parameters
Iron = {name : 'Fe', f1 : 5.26e-9, f2 : 1.24, f3 : -3.21}
Silicium = {name : 'Si',z : 1, x : 22.5, f1 : 3.77e-8, f2 : 1.48, f3 : -3.10}
Carbon = {name : 'C', z : 1, x : 30}
Oxygen = {name : 'O', z : 0.5, x : 23.2}
Sulfur = {name : 'S', z : 0.5, x : 19.4}


def rho(f,f1,f2,f3) :
    return f1*(f2 - f)**f3
    
def integrand_eq1(z):
    return z**5/(np.exp(z)-1)/(1-np.exp(-z))

def resistivity(f, T, El, theta0=417, gamma=1.52, Ir=Iron) :
    """ Ideal resistivity from the Bloch GRuneisen formula as in Gomi 2013 (equation (1))
    f: V/V0
    T : absolute temperature
    El: element (dictionnary with z valence charge and x atomic fraction)
    theta0
    gamma
    """

    theta = theta0*np.exp(-gamma *np.log(f))  #Debye temperature
    tfac = (T/theta)**5 *quad(integrand_eq1, 0, theta/T)
    
    imp = El['z']**2*El['x']

    idSi = rho(f, El['f1'], El['f2'], El['f3'])
    idFe = rho(f, Ir['f1'], Ir['f2'], Ir['f3']) * tfac

    return idFe + idSi*imp  #ideal resistivity


def k(id, T, f, L=2.44e-8):
    """ real resistivity, adding ideal resistivity and saturation resistivity
    saturation resistivity at 1bar is 1.68e-6 \Ohm m
    """
    sat = 1.68e-6 * f**(1./3.)
    return L*T*(1./id+ 1./sat)



