import numpy as np
from matplotlib import pyplot

sci_format = ":.3e"
np.set_printoptions(formatter={'float':sci_format})

# tirage aléatoire pour une loi normale de moyenne X[0] et d'incertitude-type X[1]
def tirage_normal(X):
  return X[0] + X[1]*np.random.normal()

# tirage aléatoire pour une répartition uniforme entre X[0] et X[1]
def tirage_uniforme(X):
  return X[0] + (X[1]-X[0])*np.random.random_sample()

# Accélération de la pesanteur
g = 9.80665 #en m/s, valeur normale de la CGPM

# Hauteur de chute
h1 = 3 # en m
h2 = 2 # en m
h0 = h1-h2
Deltah = 5e-3 # en m, demi-largeur de la lecture de h

# Vitesse atteinte
v0 = np.sqrt(2* g * h0) # en m/s
Deltav = 1e-2 # en m/s, incertitude-type sur la mesure de v

# Mesures
N = 10000 # nombre de points de mesure
g_calculee = np.array([])

h1min = h1-Deltah
h1max = h1+Deltah
h2min = h2-Deltah
h2max = h2+Deltah

MesuresV = np.array([tirage_normal([v0,Deltav]) for i in range(N)])
MesuresH1 = np.array([tirage_uniforme([h1min,h1max]) for i in range(N)])
MesuresH2 = np.array([tirage_uniforme([h2min,h2max]) for i in range(N)])

# Valeurs de $h = h1-h2
MesuresH = MesuresH1 - MesuresH2
# On pourrait faire la même chose sans numpy array avec:
# MesuresH = [MesuresH1[i] - MesuresH2[i] for  i in range(N)]
# ou
# for i in range(N)
#    MesuresH[i] = MesuresH1[i] - MesuresH2[i]

# Valeurs calculées de g = v^2/(2h)
g_calculees = MesuresV**2/(2*MesuresH)
g_moyenne = np.average(g_calculees)
Stdevg = np.std(g_calculees)
print('moyenne des valeurs calculées de g:',"{:.3e}".format(g_moyenne))
print('incertitude-type sur les valeurs calculées de g:',"{:.3e}".format(Stdevg))

pyplot.hist(g_calculees, bins = 50, color = 'blue', edgecolor = 'black')
pyplot.xlabel('g (m/s^2)')
pyplot.ylabel('effectif')
pyplot.title('Pour '+str(N)+ ' iterations')
# pyplot.show()

StdevRelatg = Stdevg/g_moyenne
print('écart-type relatif sur g:', "{:.3e}".format(StdevRelatg))

print('incertitude-type sur v:',Deltav)
DeltaRelatv = Deltav/v0
print('incertitude-type relative sur v:',DeltaRelatv)
Stdevh = Deltah/np.sqrt(3) # attention au facteur $\sqrt{3}$
StdevRelath1 = Stdevh/h1
StdevRelath2 = Stdevh/h2
print('incertitudes-types sur h1,h2:',Stdevh)
print('incertitudes-types relatives sur h1,h2:',StdevRelath1,StdevRelath2)

Stdevh0 = np.std(MesuresH)
print('écart-type des mesures de h1-h2' ,Stdevh0)
Deltah0 = Deltah*np.sqrt(2)/np.sqrt(3) # attention au facteur $\sqrt{3}$
print('incertitude-type composée sur h1-h2' , Deltah0)

StdevRelath0 = Stdevh0/h0
incertitude_relative_composee = np.sqrt(4*DeltaRelatv**2 + StdevRelath0**2)
print('incertitude-type relative composée sur g:',incertitude_relative_composee)
print('écart-type relatif sur g:',StdevRelatg)

Np = 100 # nombre de points de mesure
gp_calculee = np.array([])

MesuresVp = np.array([tirage_normal([v0,Deltav]) for i in range(Np)])
MesuresH1p = np.array([tirage_uniforme([h1min,h1max]) for i in range(Np)])
MesuresH2p = np.array([tirage_uniforme([h2min,h2max]) for i in range(Np)])
MesuresHp = MesuresH1p - MesuresH2p

g_calculeesp = MesuresVp**2/(2*MesuresHp)
g_moyennep = np.average(g_calculeesp)
Stdevgp = np.std(g_calculeesp)
print('moyenne des valeurs calculées de g (N=100)',g_moyennep)
print('incertitude-type sur les valeurs calculées de g (N=100):',Stdevgp)

pyplot.hist(g_calculeesp, bins = 50, color = 'blue', edgecolor = 'black')
pyplot.xlabel('g (m/s^2)')
pyplot.ylabel('effectif')
pyplot.title('Pour '+str(Np)+ ' iterations')
# pyplot.show()

m = 100
Nm = int(np.floor(N/m))

MesuresVm = np.zeros(Nm)
MesuresH1m = np.zeros(Nm)
MesuresH2m = np.zeros(Nm)
MesuresHm = np.zeros(Nm)

for i in range(Nm):
    for j in range(m):
        MesuresVm[i] += MesuresV[m*i+j]
        MesuresH1m[i] += MesuresH1[m*i+j]
        MesuresH2m[i] += MesuresH2[m*i+j]

MesuresHm = (MesuresH1m - MesuresH2m)/m
MesuresVm = MesuresVm/m

g_calculeesm = MesuresVm**2/(2*MesuresHm)
print('moyenne des moyennes de $m$ mesures', np.average(g_calculeesm))
print('écart-type des moyennes de $m$ mesures', "{:.3e}".format(np.std(g_calculeesm)))
print('écart-type des mesures brutes', "{:.3e}".format(np.std(g_calculees)))
