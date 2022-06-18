import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

figHeaviside,axHeaviside = plt.subplots()

x = np.linspace(-5,5,200)
y1 = np.piecewise(x,[x<2,x>=2],[2,3])

def Heaviside(x,params):
    x0, yavant, yapres = params
    return np.piecewise(x,[x<x0,x>=x0],[yavant,yapres])

params = [1,-4,6] 
y2 = Heaviside(x,params)

axHeaviside.set_xlabel('abscisse')
axHeaviside.set_ylabel('ordonnée')
axHeaviside.set_title('titre')
axHeaviside.plot(x,y1,'b',label='Heaviside')
axHeaviside.plot(x,y2,'r',label='Heaviside2')
axHeaviside.legend(loc='best',shadow=True)

figHeaviside.show()

NombrePoints = 2000

#  Paramètres
E = 5 # en V
R = 1e3 # en ohm
C = 1e-6 # en F

tau = R*C # en s
print(f'constante de temps: tau= {1e3*tau} ms\n valeur asymptotique: uCinf = {E} V')

tmin = -2 # en unités de tau
tmax = 6
t = np.linspace(tmin,tmax,NombrePoints)
Deltat = (tmax-tmin)/(NombrePoints -1)

def dxoverdt(t,x,params):
    t0,tau,xinfty = params
    return np.piecewise(t,[t<t0,t>=t0],[-x/tau,(xinfty-x)/tau])

x = np.zeros(NombrePoints)

params = [0,1,3]

for i in range(NombrePoints-1):
    x[i+1] = x[i]+ Deltat*dxoverdt(t[i],x[i],params)

figequadiff,axequadiff = plt.subplots()

axequadiff.set_xlabel('t/tau')
axequadiff.set_ylabel('u/E')
axequadiff.set_title('Réponse d\'un dipôle RC à un échelon de tension')
axequadiff.plot(t,x,'b',label='uc/E')
axequadiff.legend(loc='best',shadow=True)

figequadiff.show()

from scipy.integrate import solve_ivp

# paramètres de l'échelon
t0 = 0
xavant = 2
xapres = 5

def dxoverdt(t,x,t0,xavant,xapres):
    return np.piecewise(t,[t<t0,t>=t0],[xavant-x,xapres-x])

# intervalle de temps
tmin = -2 # en unités de tau
tmax = 5

# évènements à rechercher: ici passage par x=4
def passage(t,x,t0,xavant,xapres):
    return x[0]-4

# Condition initiale (à t = tmin)
CI = [0]

# Résolution numérique de l'équation différentielle
solution = solve_ivp(dxoverdt,[tmin,tmax],CI,args=(t0,xavant,xapres),events=passage,dense_output=True,max_step=.1)

print(f'passage à {solution.y_events} en t={solution.t_events}')

figSolveIVP,axSolveIVP = plt.subplots()

axSolveIVP.set_xlabel('t/tau')
axSolveIVP.set_ylabel('u/E')
axSolveIVP.set_title('Résolution par SolveIVP')
axSolveIVP.plot(solution.t,solution.y[0],'b+',label='Valeurs évaluées')

NombrePoints = 200
instants = np.linspace(tmin,tmax,NombrePoints)
valeurs = solution.sol(instants)
axSolveIVP.plot(instants,valeurs.T,'r',label='Fonction interpolée')
axSolveIVP.legend(loc='best',shadow=True)

figSolveIVP.show()

NombrePoints = 2000

#  Paramètres
E = 10 # en V
R = 500 # en ohm
L = 2.9e-2 # en H

uR = E # en V
tau = L/R # en s
print(f'constante de temps: tau= {1e3*tau} ms\n valeur asymptotique: uR = {uR} V')

tmin = -.2 # en unités de tau
tmax = 6
t = np.linspace(tmin,tmax,NombrePoints)
Deltat = (tmax-tmin)/(NombrePoints -1)

def dxoverdt(t,x,params):
    t0,tau,xinfty = params
    return np.piecewise(t,[t<t0,t>=t0],[-x/tau,(xinfty-x)/tau])

x = np.zeros(NombrePoints)

params = [0,1,E]

for i in range(NombrePoints-1):
    x[i+1] = x[i]+ Deltat*dxoverdt(t[i],x[i],params)

figequadiff,axequadiff = plt.subplots()

axequadiff.set_xlabel('t/tau')
axequadiff.set_ylabel('uR (V)')
axequadiff.set_title('Réponse d\'un dipôle RL à un échelon de tension')
axequadiff.plot(t,x,'b',label='uR')
axequadiff.legend(loc='best',shadow=True)

figequadiff.show()

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
figzoom,axezoom = plt.subplots()
axezoom.xaxis.set_minor_locator(MultipleLocator(0.01))
axezoom.yaxis.set_minor_locator(MultipleLocator(0.005))

axezoom.set_xlabel('t/tau')
axezoom.set_ylabel('uR (V)')
axezoom.set_title('Réponse d\'un dipôle RL à un échelon de tension')
axezoom.plot(t,x,'b',label='uR')
axezoom.legend(loc='best',shadow=True)
axezoom.set_xlim(.9,.93)
axezoom.set_ylim(5.8,6.2)
axezoom.grid(which='both')

figzoom.show()

def dxoverdtCommutation(t,x,params):
    periode, tau, ybas, yhaut = params
    demiperiode = periode/2
    parite = (np.floor(t//demiperiode))%2 # vaut 0 ou 1 selon la 1/2 période dans laquelle on est
    return (np.piecewise(t,[parite != 0, parite == 0 ],[ybas,yhaut]) -x)/tau

iLinf = 1000*E/R # en mA
tau = 1000*L/R # en ms
tmin = 0
periode = 1/5 # en ms
NombrePeriodes = 10
tmax = NombrePeriodes * periode
t = np.linspace(tmin,tmax,NombrePoints)
Deltat = (tmax-tmin)/(NombrePoints -1)
params = [periode,tau,-iLinf,iLinf] # t ms
x = np.zeros(NombrePoints)

for i in range(NombrePoints-1):
    x[i+1] = x[i]+ Deltat*dxoverdtCommutation(t[i],x[i],params)

figCommutation,axCommutation = plt.subplots()

axCommutation.set_xlabel('t (ms)')
axCommutation.set_ylabel('iL (mA)')
axCommutation.set_title('Réponse d\'un dipôle RL à un créneau de tension')
axCommutation.plot(t,x,'b',label='iL')
axCommutation.legend(loc='best',shadow=True)

figCommutation.show()

NombrePointsDemiPeriode = int((periode/(2*Deltat)))+1

iLmoyenneA = np.zeros(NombrePeriodes) # 1/2 périodes impaires
iLminA = np.zeros(NombrePeriodes)
iLmaxA = np.zeros(NombrePeriodes)

for i in range(NombrePeriodes):
    iLmoyenneA[i] = np.average(x[2*i*NombrePointsDemiPeriode:(2*i+1)*NombrePointsDemiPeriode])
    iLminA[i] = np.min(x[2*i*NombrePointsDemiPeriode:(2*i+1)*NombrePointsDemiPeriode])
    iLmaxA[i] = np.max(x[2*i*NombrePointsDemiPeriode:(2*i+1)*NombrePointsDemiPeriode])

print('1/2 périodes impaires')
print(f'valeurs moyennes de iL (en mA): {iLmoyenneA}\nvaleurs max: {iLmaxA}\nvaleurs min {iLminA}')


iLmoyenneB = np.zeros(NombrePeriodes) # 1/2 périodes paires
iLminB = np.zeros(NombrePeriodes)
iLmaxB = np.zeros(NombrePeriodes)

for i in range(NombrePeriodes):
    iLmoyenneB[i] = np.average(x[(2*i+1)*NombrePointsDemiPeriode:(2*i+2)*NombrePointsDemiPeriode])
    iLminB[i] = np.min(x[(2*i+1)*NombrePointsDemiPeriode:(2*i+2)*NombrePointsDemiPeriode])
    iLmaxB[i] = np.max(x[(2*i+1)*NombrePointsDemiPeriode:(2*i+2)*NombrePointsDemiPeriode])

print('1/2 périodes paires')
print(f'valeurs moyennes de iL (en mA): {iLmoyenneB}\nvaleurs max: {iLmaxB}\nvaleurs min {iLminB}')

iLmoyenneC = np.zeros(NombrePeriodes)
iLminC = np.zeros(NombrePeriodes)
iLmaxC = np.zeros(NombrePeriodes)
for i in range(NombrePeriodes):
    xTemp = [x[j] for j in range(len(x)) if t[j] >= i*periode and t[j] < (i+1/2) * periode]
    iLmoyenneC[i] = np.average(xTemp)
    iLminC[i] = np.min(xTemp)
    iLmaxC[i] = np.max(xTemp)

print('1/2 periodes impaires')
print(f'valeurs moyennes de iL (en mA): {iLmoyenneC}\nvaleurs max: {iLmaxC}\nvaleurs min {iLminC}')

def dxoverdtTriangle(t,x,params):
  periode, tau, ybas, yhaut = params
  demiperiode = periode/2
  reste = t % demiperiode
  pente = (yhaut-ybas)/demiperiode
  parite = (t//demiperiode)%2 # vaut 0 ou 1 selon la 1/2 période dans laquelle on est
  return  (np.piecewise(t,[(t//demiperiode)%2 != 0, (t//demiperiode)%2 == 0 ],[lambda t: ybas+pente*(t % demiperiode), lambda t: yhaut-pente*(t % demiperiode)])-x)/tau

x = np.zeros(NombrePoints)

periode = 1/5 # en ms
tau = 1000*L/R # en ms
params = [periode,tau,-E/L,E/L] # E/L en A/s soit en mA/ms

for i in range(NombrePoints-1):
    x[i+1] = x[i]+ Deltat*dxoverdtTriangle(t[i],x[i],params)

figTriangle,axTriangle = plt.subplots()

axTriangle.set_xlabel('t (ms)')
axTriangle.set_ylabel('iL (mA)')
axTriangle.set_title('Réponse d\'un dipôle RL à une fonction triangle')
axTriangle.plot(t,x,'b',label='iL')
axTriangle.legend(loc='best',shadow=True)

figTriangle.show()
