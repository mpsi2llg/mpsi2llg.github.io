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

params = [0,1,E]

for i in range(NombrePoints-1):
    x[i+1] = x[i]+ Deltat*dxoverdt(t[i],x[i],params)

figequadiff,axequadiff = plt.subplots()

axequadiff.set_xlabel('t/tau')
axequadiff.set_ylabel('u/E')
axequadiff.set_title('Réponse d\'un dipôle RC à un échelon de tension')
axequadiff.plot(t,x,'b',label='uc')
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
