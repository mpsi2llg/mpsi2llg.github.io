import numpy as np
from scipy.optimize import fsolve

def equation(x):
    return x**5 - 1e-2*np.sqrt(2-x)

root = fsolve(equation,.2)
print(f'La solution est : {root}')
