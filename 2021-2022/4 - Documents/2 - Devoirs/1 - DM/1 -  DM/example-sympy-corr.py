from sympy import linsolve, pprint, symbols

i, i1, i2, i3, u1, u2 = symbols('i i1 i2 i3 u1 u2')  # Inconnues
E, R, Rf = symbols('E R Rf')  # Données
systeme = [i1 + i2 + i3 -i,
           u1 + 2*R*i -E,
           u2 - R *i2,
           u2 -(2*R + Rf)*i3,
           u1 -R*i1,
           u1-u2 - 2*R*(i2+i3)
           ]

solution = linsolve(systeme,i, i1,i2,i3,u1,u2)
pprint(solution)
