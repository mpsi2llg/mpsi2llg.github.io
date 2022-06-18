import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def systdiff(tau,y):
    theta,thetaprime = y
    # d theta/d t = thetaprime
    # d thetaprime / dt = - sin(theta)
    return [thetaprime, - (2*np.pi)**2*np.sin(theta)]

def passage_origine(tau,y):
    theta,thetaprime=y
    return theta
passage_origine.terminal = False #pour poursuivre l'intégration
passage_origine.direction = 1 #pour ne compter que les passages avec theta croissant

longueur = .4 #m
g0 = 9.8 #m/s^2
omega0 = np.sqrt(g0/longueur) #rad/s
T0 = 2*np.pi/omega0

tau_min = 0
tau_max = 5 #périodes T0

theta0 = -np.pi/2 #angle initial (rad)
v0 =  0 #vitesse (m/s)
thetaprime0 = v0/(longueur*T0) # (rad)
CI = [theta0,thetaprime0]

pendule = solve_ivp(systdiff,[tau_min,tau_max], CI,max_step= T0/50,events= passage_origine)
angles = pendule.y[0] #en rad
anglesDeg = angles*180/np.pi #en deg
instantsAdim = pendule.t #en unités de T0
instants = instantsAdim*T0 #en s
vitessesAngAdim = pendule.y[1] #en unités de 1/T0
vitesses = vitessesAngAdim*longueur/T0 #en m/s

f'période pour theta0 = {theta0*180/np.pi} deg: {np.mean(np.diff(pendule.t_events[0]*T0)):.2E} s' #np.diff calcule la différence des termes consécutifs de la liste

f'période des petits angle : {T0:.2E} s'

fig,(axtemp,axphase) = plt.subplots(1,2)
fig.tight_layout()
axtemp.plot(instants,angles)
axphase.plot(angles,vitesses)
axtemp.set_xlabel(r"$t$ (s)")
axtemp.set_ylabel(r"$\theta$ (deg)")
axphase.set_xlabel(r"$\theta$ (deg)")
axphase.set_ylabel(r"$v$ (m/s)")
axphase.yaxis.set_label_position("right")
fig.show()

def systdiffDM(tau,y,omega2g,omega2k):
    alpha,alphaprime = y
    # d alpha/d t = alphaprime
    # d alphaprime / dt = - (g/R) cos(alpha) + (k/m) (sin(alpha) - sin(alpha/2))
    return [alphaprime, - omega2g * np.cos(alpha) + omega2k * ( np.sin (alpha) - np.sin(alpha/2))]
def passage_equilibre(tau,y,omega2g,omega2k):
      alpha,alphaprime=y
      return alpha-3*np.pi/4
passage_equilibre.terminal = False #pour poursuivre l'intégration
passage_equilibre.direction = 1 #pour ne compter que les passages avec theta croissant

g = 9.8 #m/s^2
R = .1 #m
m = 5e-2 #kg
k = 15.6 #N/m
omegag = np.sqrt(g/R) # rad/s
omegak = np.sqrt(k/m) # rad/s
omega2g = omegag**2
omega2k = omegak**2

omegaDM = np.sqrt((omega2g - omega2k) * np.cos(3*np.pi/4) + (omega2k/2)* np.cos(3*np.pi/8))
TDM = 2*np.pi/omegaDM

alpha03a = 3*np.pi/4 - np.pi/10
alphaprime03a = 0

tmin3a,tmax3a = 0,5*TDM
CI3a = [alpha03a, alphaprime03a]

mouvement3a = solve_ivp(systdiffDM,[tmin3a,tmax3a],CI3a,max_step= TDM/50, args=[omega2g,omega2k],events=passage_equilibre)
angles3a = mouvement3a.y[0] #en rad
anglesDeg3a = angles3a*180/np.pi #en deg
instants3a= mouvement3a.t #en s
vitesses3aAng = mouvement3a.y[1] #en rad/s
vitesses3a = vitesses3aAng*R #en m/s

fig3a,(ax3atemp,ax3aphase) = plt.subplots(1,2)
fig3a.tight_layout()
ax3atemp.plot(instants3a,anglesDeg3a)
ax3aphase.plot(anglesDeg3a,vitesses3a)
ax3atemp.set_xlabel(r"$t$ (s)")
ax3atemp.set_ylabel(r"$\alpha$ (deg)")
ax3aphase.set_xlabel(r"$\alpha$ (deg)")
ax3aphase.set_ylabel(r"$v$ (m/s)")
ax3aphase.yaxis.set_label_position("right")
fig3a.show()

f'période des petites oscillations: {TDM:.2E} s' #np.diff calcule la différence des termes consécutifs de la liste

f'période pour alpha0 = {alpha03a*180/np.pi} deg: {np.mean(np.diff(mouvement3a.t_events[0])):.2E} s' #np.diff calcule la différence des termes consécutifs de la liste

alpha03b = np.pi
alphaprime03b = 0

tmin3b,tmax3b = 0,5*TDM
CI3b = [alpha03b, alphaprime03b]

mouvement3b = solve_ivp(systdiffDM,[tmin3b,tmax3b],CI3b,max_step= TDM/50, args=[omega2g,omega2k],events=passage_equilibre)
angles3b = mouvement3b.y[0] #en rad
anglesDeg3b = angles3b*180/np.pi #en deg
instants3b= mouvement3b.t #en s
vitesses3bAng = mouvement3b.y[1] #en rad/s
vitesses3b = vitesses3bAng*R #en m/s

fig3b,(ax3btemp,ax3bphase) = plt.subplots(1,2)
fig3b.tight_layout()
ax3btemp.plot(instants3b,anglesDeg3b)
ax3bphase.plot(anglesDeg3b,vitesses3b)
ax3btemp.set_xlabel(r"$t$ (s)")
ax3btemp.set_ylabel(r"$\alpha$ (deg)")
ax3bphase.set_xlabel(r"$\alpha$ (deg)")
ax3bphase.set_ylabel(r"$v$ (m/s)")
ax3bphase.yaxis.set_label_position("right")
fig3b.show()

f'période pour alpha0 = {alpha03b*180/np.pi} deg: {np.mean(np.diff(mouvement3b.t_events[0])):.2E} s' #np.diff calcule la différence des termes consécutifs de la liste

f'écart relatif de la période pour alpha0 = {alpha03b*180/np.pi} deg: {100*(np.mean(np.diff(mouvement3b.t_events[0]))/TDM -1):.2E} %'

alpha03c = np.pi/4
alphaprime03c = 0

tmin3c,tmax3c = 0,5*TDM
CI3c = [alpha03c, alphaprime03c]

def passage_A(tau,y,omega2g,omega2k):
    alpha,alphaprime=y
    return alpha-np.pi
passage_A.terminal = True #pour poursuivre l'intégration

mouvement3c = solve_ivp(systdiffDM,[tmin3c,tmax3c],CI3c,max_step= TDM/50, args=[omega2g,omega2k],events=passage_A)
angles3c = mouvement3c.y[0] #en rad
anglesDeg3c = angles3c*180/np.pi #en deg
instants3c= mouvement3c.t #en s
vitesses3cAng = mouvement3c.y[1] #en rad/s
vitesses3c = vitesses3cAng*R #en m/s

fig3c,(ax3ctemp,ax3cphase) = plt.subplots(1,2)
fig3c.tight_layout()
ax3ctemp.plot(instants3c,anglesDeg3c)
ax3cphase.plot(anglesDeg3c,vitesses3c)
ax3ctemp.set_xlabel(r"$t$ (s)")
ax3ctemp.set_ylabel(r"$\alpha$ (deg)")
ax3cphase.set_xlabel(r"$\alpha$ (deg)")
ax3cphase.set_ylabel(r"$v$ (m/s)")
ax3cphase.yaxis.set_label_position("right")
fig3c.show()

f'vitesse quand il parvient en B: {mouvement3c.y_events[0][0][0]:.2E} m/s'

fig3d,ax3dphase = plt.subplots()
fig3d.tight_layout()
ax3dphase.plot(anglesDeg3a,vitesses3a)
ax3dphase.plot(anglesDeg3c,vitesses3c)
ax3dphase.plot(anglesDeg3b,vitesses3b)
ax3dphase.set_xlabel(r"$\alpha$ (deg)")
ax3dphase.set_ylabel(r"$v$ (m/s)")
ax3dphase.yaxis.set_label_position("right")
fig3d.show()
