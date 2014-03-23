import matplotlib.pyplot as plt
import sympy
import numpy as np

roots = ((0, 0), (2, 0), (0, 2), (2, 2))
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
            e_str = '%d' if degree < 10 else '_%d'
            symbols.append('c%s' % ''.join(e_str % e for e in exps))
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

fp = sympy.lambdify(xs, p, modules='numpy')
print("Polynomium:")
print(sympy.Eq(p, 0))
for xi in roots:
    print("Check root %s:" % (xi,))
    print(fp(*xi))

x1 = min(root[0] for root in roots)
x2 = max(root[0] for root in roots)
y1 = min(root[1] for root in roots)
y2 = max(root[1] for root in roots)

w = x2 - x1
h = y2 - y1

if w < 1e-6:
    w = 1
if h < 1e-6:
    h = 1

x1 -= w/10
x2 += w/10
y1 -= h/10
y2 += h/10

fig = plt.figure()
ax = fig.add_subplot(111)
xs = np.linspace(x1, x2, 200)
ys = np.linspace(y1, y2, 200)
xy = np.meshgrid(xs, ys)
ax.pcolormesh(xy[0], xy[1], np.minimum(1, np.maximum(-1, fp(xy[0], xy[1]))))
rootx, rooty = zip(*roots)
ax.plot(rootx, rooty, 'o')

plt.show()
