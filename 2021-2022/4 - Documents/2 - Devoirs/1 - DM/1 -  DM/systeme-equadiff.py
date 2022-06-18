import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy.ma as ma

def systdiff(u,tau):
    theta,thetaprime = u
    # d theta/d t = thetaprime
    # d thetaprime / dt = - sin(theta)
    return [thetaprime, - (2*np.pi)**2*np.sin(theta)]

longueur = .4 #m
g0 = 9.8 #m/s^2
omega0 = np.sqrt(g0/longueur) #rad/s
T0 = 2*np.pi/omega0

tau_min = 0
tau_max = 5 #périodes T0
NombrePoints = 2000
tau = np.linspace(tau_min,tau_max,NombrePoints)
mask=ma.masked_greater(tau,1).mask #pour ne conserver que l'intervalle tau = 0:1, soit t = 0:T0
instants = tau*T0
instantsMasked = instants[~mask]

theta0 = np.pi/2 #angle initial (rad)
v0 =  2 #vitesse (m/s)
thetaprime0 = v0/(longueur*T0) # (rad)
CI = [theta0,thetaprime0]

sol = odeint(systdiff,CI,tau)
angles = sol[:,0]
vitessesAngAdim = sol[:,1] #en unités de 1/T_0
vitessesAng = vitessesAngAdim/T0 #en unités de 1/T_0
vitesses = 100*vitessesAngAdim*longueur/T0 # en cm/s
anglesMasked = angles[~mask]
vitessesAngMasked = vitessesAng[~mask]
vitessesMasked = vitesses[~mask]

figtemporel,(axtempAdim,axtempDim) = plt.subplots(1,2) #pour avoir deux figures côte à côte
figtemporel.tight_layout()
axtempAdim.plot(tau,angles,label='Angle')
axtempAdim.set_xlabel(r"$t/T_0$")
axtempAdim.set_ylabel(r"$\tau(rad)$" )
axtempAdim.legend(loc='best',shadow=True)

axtempDim.plot(instants,angles,label='Angle')
axtempDim.set_xlabel(r"$t(s)$")
axtempDim.set_ylabel(r"$\tau(rad)$" )
axtempDim.yaxis.set_label_position("right")
axtempDim.legend(loc='best',shadow=True)

figtemporel.show()

figphase,axphase = plt.subplots()
figphase.tight_layout()
axphase.plot(anglesMasked,vitessesMasked)
axphase.set_xlabel(r"$\theta$ (rad)")
axphase.set_ylabel(r"$v_\theta$(m/s)" )

figphase.show()
