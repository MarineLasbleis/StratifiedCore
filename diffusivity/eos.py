#!/usr/bin/env python
## Time-stamp: "2015-12-03 16:57:05 marine"

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad




def Volume(P, T, param):
    """ Calculation of the typical volume
    
    P : pressure, T : temperature
    param : parameters for the Burch Murnaghan EOS
    
     """
    rhodP = param['rho0'] * (param['KPrim'] * P/param['K'] + 1.)**(1./param['KPrim'])
    rhodT = np.exp( - (param['alpha'] * np.exp(-param['alphaPrim'] *P / param['K']))
                    * (T-param['T0']) )

    return rhodP*rhodT




def fval(P,T, param):
    """ V0/V """
    return Volume(0., 300., param)/Volume(P, T, param)



def integrande_eq(p, t, param):
    """ used in the function Tadibatic"""
    alpha = param['alpha'] * np.exp( -param['alphaPrim']/param['K']*p)
    return alpha/Volume(p, t, param)

def Tadiabatic(P, t0, p0, param):
    """ Compute the adiabatic profile from (T0,P0) to P"""
   
    intexp = quad(integrande_eq, p0, P, args=(t0,param))
    intexp = intexp[0]
    return t0*np.exp( intexp /param['Cp'] )
    





### Test
if __name__ == '__main__':


    MURNAGHAN = {'T0': 1812., 'rho0': 7010., 'K': 130.*1e9, 'KPrim':4.,
                        'alphaPrim': 130.*np.log(2.)/(360.-135.) ,
                        'alpha': 1.e-5/np.exp(-0.400485 *135./130.),   #1.e-5*np.exp(130.*np.log(2.)/(360.-135.)*135./130),#   1.e-5/np.exp( np.log(2)/(360.-135.)*135.),
                        'Cp': 800., 'cond': 140.}


    
    P = np.linspace(135.,330.,100)*1.e9
    T = []

    P0 = 135e9
    T0 = 3750.
    
    for p in P:
        #print Tadiabatic(p, T0, P0, MURNAGHAN)
        T = np.append(T, Tadiabatic(p, T0, P0, MURNAGHAN))


    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(P/1.e9,T)
    ax.annotate("("+str(P0/1e9)+","+str(T0)+")", xy=(P0/1e9, T0), xycoords='data',
                xytext=(-50, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
    ax.annotate("("+str(P[-1]/1e9)+","+str(T[-1])+")", xy=(P[-1]/1e9, T[-1]), xycoords='data',
                xytext=(-50, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
    ax.set_xlabel("Pressure (GPa)")
    ax.set_xlabel("Temperature (K)")
    ax.set_title("Adiabatic profile")
    
    plt.show()
    
    

