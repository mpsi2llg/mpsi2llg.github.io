from sympy import linsolve, pprint, symbols

x, y, z = symbols('x y z')  # Inconnues
a1, a2, a3 = symbols('a1 a2 a3')  # Données
systeme = [x + y +z -a1,
           2* x + y -z -a2,
           x - 2*y-a3]

solution = linsolve(systeme,x,y,z)
pprint(solution)
