import numpy as np
from scipy.optimize import fsolve

def equation(x):
    return x**5 - 1e-2*np.sqrt(2-x)

root = fsolve(equation,.2)
print(f'La solution est : {root}')

c = 1e-2 # en mol.L
K1 = 10**(8.3)
K3 = 10**(15.2)

def equation3(x):
    return x**5 - np.sqrt(K3)/(c*K1)*(1-x)**(3/2)

root = fsolve(equation3,.87)
print(f'Pour c = {c}, la solution est du 3a est: tau={root}')

c = 1e-1 # en mol.L
root = fsolve(equation3,.87)
print(f'Pour c = {c}, la solution est du 3a est: tau={root}')
