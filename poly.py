import sympy

roots = ((-1, -1), (0, 0), (1, 2))
degree = 2
dims = 2

if dims <= 3:
    names = ('x', 'y', 'z')
    xs = sympy.symbols(names[:dims])
else:
    xs = sympy.symbols(tuple('x%d' % n for n in range(1, dims+1)))

def sum_to(n, degree):
    if n == 1:
        yield (degree,)
        return
    for d in range(degree+1):
        for p in sum_to(n-1, degree-d):
            yield (d,) + p

def monomials(xs, degree):
    for d in range(degree+1):
        for exps in sum_to(len(xs), d):
            p = 1
            for x, e in zip(xs, exps):
                p *= x ** e
            yield p

def coefficients(n, degree):
    symbols = []
    for d in range(degree+1):
        for exps in sum_to(n, degree):
            symbols.append('c%s' % ''.join('_%d' % e for e in exps))
    return tuple(sympy.symbols(symbols))

a = tuple(tuple(monomials(x, degree)) for x in roots)
b = tuple(0 for each in roots)
cs = tuple(coefficients(dims, degree))

eqns = tuple(sympy.Eq(sum(ci * ai for ci, ai in zip(cs, row)), bi) for row, bi in zip(a, b))
print("Linear system:")
print(eqns)
sol = sympy.solve(eqns, cs)
print("Solution(s):")
print(sol)

p = sum(sol.get(ci, ci) * m for ci, m in zip(cs, monomials(xs, degree)))
p = p.subs((ci, 1) for ci in cs)
print("Polynomium:")
print(sympy.Eq(p, 0))
for xi in roots:
    print("Check root %s:" % (xi,))
    print(p.subs(zip(xs, xi)))
