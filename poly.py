import sympy

roots = (-1, 1, 4, 6)
degree = 4

a = tuple(tuple(x ** d for d in range(degree+1)) for x in roots)
b = tuple(0 for each in range(degree+1))
cs = tuple(sympy.symbols('c%d' % d) for d in range(degree+1))

eqns = tuple(sympy.Eq(sum(ci * ai for ci, ai in zip(cs, row)), bi) for row, bi in zip(a, b))
print(eqns)
sol = sympy.solve(eqns, cs)

x = sympy.symbols('x')
p = sum(sol.get(ci, ci) * x ** d for ci, d in zip(cs, range(degree+1)))
print(sympy.Eq(p, 0))
for xi in roots:
    print("Check root %d:" % xi)
    print(p.subs(x, xi))
