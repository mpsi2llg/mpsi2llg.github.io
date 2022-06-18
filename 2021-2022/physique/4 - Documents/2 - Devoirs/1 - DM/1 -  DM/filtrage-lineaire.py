import numpy as np
from scipy import signal
from scipy import fftpack
import matplotlib.pyplot as plt
import numpy.ma as ma

pi = np.pi
# Définition du filtre
f_c = 2e2              # fréquence de coupure (en Hz)
omega_c = 2.0*pi*f_c   # pulsation de coupure (en rad/s)
H_0 = 2                # gain en bande passante   
# Coefficients du dénominateur rangés par degrés décroissants
a = np.array( [1/omega_c,1.0] )
# Coefficients du numérateur
b = np.array( [H_0] )
H_coeffs = (b,a)

# Plage de fréquence à étudier
frequences = np.logspace(1 , 5 , 400)
[ w, H ] = signal.freqs(b,a,2*pi*frequences)

figBode1,(axBodeGain1,axBodePhase1) = plt.subplots(1,2) #pour avoir deux figures côte à côte
figBode1.tight_layout()
# gain en décibel, np.absolute calcule le module
G = 20.0 * np.log10(np.absolute(H))
# phase en degrés, np.angle calcule l'argument en radian, converti par rad2deg
phase = np.rad2deg(np.angle(H))
axBodeGain1.semilogx()
axBodePhase1.semilogx()
axBodeGain1.plot(frequences,G,label='Gain')
axBodeGain1.set_xlabel(r"$f(Hz)$")
axBodeGain1.set_ylabel(r"$G_{dB}$" )
axBodePhase1.plot(frequences,phase,label='Phase')
axBodePhase1.set_xlabel(r"$f(Hz)$")
axBodePhase1.set_ylabel(r"$phi(°)$")
axBodePhase1.yaxis.set_label_position("right")
axBodePhase1.legend(loc='best',shadow=True)
axBodeGain1.legend(loc='best',shadow=True)
figBode1.show()

HarmoniqueMax = 9
fCre = 100 #Hz
Delta_t = 1/(20*HarmoniqueMax*fCre)
echantillon_t = np.arange(0,20/fCre,Delta_t) # pour disposer de 20 périodes pour l'échantillonage
EntreeCreneau = signal.square(echantillon_t*2*pi*fCre) # signal.square a pour période 2 pi
figCreneau,axCreneau = plt.subplots()
axCreneau.set_xlim(0,3/fCre)
axCreneau.plot(echantillon_t,EntreeCreneau)
axCreneau.set_xlabel(r"$t(s)$")
axCreneau.set_ylabel(r"U(V)" )
figCreneau.show()

FrequencesFFTCreneau =  fftpack.fftfreq(len(EntreeCreneau) , Delta_t)
FFTEntreeCreneau = fftpack.fft(EntreeCreneau)
CutoffGdB=40
FFTEntreeCreneauGdB=20*np.log10(np.absolute(FFTEntreeCreneau))
mask=ma.masked_less(FFTEntreeCreneauGdB,np.max(FFTEntreeCreneauGdB)-CutoffGdB).mask
FrequencesFFTCreneauMasked=FrequencesFFTCreneau[~mask]
FFTEntreeCreneauMasked=FFTEntreeCreneau[~mask]
FFTEntreeCreneauGdBMasked=FFTEntreeCreneauGdB[~mask]
figFFTCreneau,(axFFTEntreeCreneau,axFFTSortieCreneau) = plt.subplots(1,2) #pour avoir deux figures côte à côte
figFFTCreneau.tight_layout()
axFFTEntreeCreneau.plot(FrequencesFFTCreneauMasked,np.absolute(FFTEntreeCreneauMasked),'o',label='Entree')
axFFTEntreeCreneau.set_xlabel(r"$f(Hz)$")
axFFTEntreeCreneau.set_ylabel(r"amplitude entrée" )
axFFTEntreeCreneau.set_xlim([0, HarmoniqueMax*fCre])
axFFTEntreeCreneau.grid(which='both')
axFFTEntreeCreneau.legend(loc='best',shadow=True)

HFFTCreneau = signal.freqs(b,a,2*pi*FrequencesFFTCreneau)[1]
FFTSortieCreneau = HFFTCreneau * FFTEntreeCreneau
FFTSortieCreneauMasked = FFTSortieCreneau[~mask]
axFFTSortieCreneau.plot(FrequencesFFTCreneauMasked,np.absolute(FFTSortieCreneauMasked),'+',label='Sortie')
axFFTSortieCreneau.set_xlabel(r"$f(Hz)$")
axFFTSortieCreneau.set_ylabel(r"amplitude sortie" )
axFFTSortieCreneau.set_xlim([0, HarmoniqueMax*fCre])
axFFTSortieCreneau.grid(which='both')
axFFTSortieCreneau.yaxis.set_label_position("right")
axFFTSortieCreneau.legend(loc='best',shadow=True)
figFFTCreneau.show()

SortieCreneau = fftpack.ifft(FFTSortieCreneau)
figCreneauES,axCreneauES = plt.subplots()
axCreneauES.plot(echantillon_t,EntreeCreneau,label='entree')
print(f'max des parties imaginaires {max(np.imag(SortieCreneau))}')
axCreneauES.plot(echantillon_t,np.real(SortieCreneau),label='sortie') 
axCreneauES.set_xlim([0, 3/fCre])
axCreneauES.set_xlabel(r't(s)')
axCreneauES.set_ylabel(r'U(V)')
axCreneauES.legend(loc='best',shadow=True)
figCreneauES.show()
