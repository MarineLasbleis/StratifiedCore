#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


ParametersMurnaghan = {'T0': 1812., 'rho0': 7010., 'K': 130.*1e9, 'KPrim':4.,
                        'alphaPrim': 130.*np.log(2.)/(360.-135.) ,
                        'alpha': 1.e-5/np.exp(-0.400485 *135./130.),   #1.e-5*np.exp(130.*np.log(2.)/(360.-135.)*135./130),#   1.e-5/np.exp( np.log(2)/(360.-135.)*135.),
                        'Cp': 800., 'cond': 140.}


def Volume(P, T, param):
    """ Calculation of the typical volume
    P : pressure, T : temperature
    param : parameters for the Burch Murnaghan EOS
     """
    rhodP = param['rho0'] * (param['KPrim'] * P/param['K'] + 1.)**(1./param['KPrim'])
    rhodT = np.exp( - (param['alpha'] * np.exp(-param['alphaPrim'] *P / param['K']))
                    * (T-param['T0']) )

    return rhodP*rhodT

def integrande_eq(p, t, param):
    alpha = param['alpha'] * np.exp( -param['alphaPrim']/param['K']*p)
    return alpha/Volume(p, t, param)    

def Tadiabatic(P, t0, p0, param):
    """ Compute the adiabatic profile from (T0,P0) to P"""
   
    intexp = quad(integrande_eq, p0, P, args=(t0,param))
    intexp = intexp[0]
    return t0*np.exp( intexp /param['Cp'] )
    
def fval(P,T, param):
    """ V0/V """
    return Volume(0., 300., param)/Volume(P, T, param)




### Test
if __name__ == '__main__':
    P = np.linspace(135.,330.,100)*1.e9
    T = []

    print 'alpha', ParametersMurnaghan['alpha'], ParametersMurnaghan['alphaPrim']

    P0 = 135e9
    T0 = 3750.
    
    for p in P:
        print Tadiabatic(p, T0, P0, ParametersMurnaghan)
        T0 = Tadiabatic(p, T0, P0, ParametersMurnaghan)
        T = np.append(T, T0)
        P0=p
        P0 = P[0]
        T0 = 3750.

    print P[0], len(P), len(T)
    plt.plot(P/1.e9,T)

    print Volume(135., 3750., ParametersMurnaghan)

    print ParametersMurnaghan
    
    plt.show()
    
    

