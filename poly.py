import matplotlib.pyplot as plt
import sympy
import numpy as np

roots = ((0, 0), (2, 1), (1, 2), (3, 3))
dims = 2

def default_dimensions(n):
    """Default names (in order) given to axes in n dimensional space."""
    if n <= 3:
        names = ('x', 'y', 'z')
        return sympy.symbols(names[:n])
    else:
        return sympy.symbols(tuple('x%d' % n for n in range(1, n+1)))

xs = default_dimensions(dims)

def sum_to(n, degree):
    """Yield all ways to sum to `degree` using `n` integers."""
    if n == 1:
        yield (degree,)
        return
    for d in range(degree+1):
        for p in sum_to(n-1, degree-d):
            yield (d,) + p

def monomials(xs, degree):
    """Yield all monomials over `xs` of degree at most `degree`."""
    for d in range(degree+1):
        for exps in sum_to(len(xs), d):
            p = 1
            for x, e in zip(xs, exps):
                p *= x ** e
            yield p

def coefficients(n, degree):
    symbols = []
    for d in range(degree+1):
        for exps in sum_to(n, d):
            e_str = '%d' if d < 10 else '_%d'
            symbols.append('c%s' % ''.join(e_str % e for e in exps))
    return tuple(sympy.symbols(symbols))

def polynomial_through(roots, degree, xs, *, general=False, verbose=False):
    a = tuple(tuple(monomials(x, degree)) for x in roots)
    b = tuple(0 for each in roots)
    cs = tuple(coefficients(dims, degree))

    eqns = tuple(sympy.Eq(sum(ci * ai for ci, ai in zip(cs, row)), bi) for row, bi in zip(a, b))
    if verbose:
        print("Linear system:")
        print(eqns)
    sol = sympy.solve(eqns, cs)
    if verbose:
        print("Solution(s):")
        print(sol)

    if () == tuple(filter(lambda n: n != 0, sol.values())):
        if verbose:
            print("No non-trivial solution")
        return None

    p = sum(sol.get(ci, ci) * m for ci, m in zip(cs, monomials(xs, degree)))
    frees = tuple(frozenset(p.free_symbols).difference(xs))

    if verbose:
        print("Free symbols:")
        print(frees)

    if not general:
        concretes = []
        for free in frees:
            concretes.append(p.subs((v, 1 if v == free else 0) for v in frees))
        degrees = map(lambda p: sympy.Poly(p, *xs).total_degree(), concretes)
        ps, ds = zip(*sorted(zip(concretes, degrees), key=lambda a: a[1]))
        if verbose:
            print(ps)
            print(ds)
        p = ps[0]

    return p

max_degree = int(1 + dims * len(roots) ** (1/dims))
min_degree = 1
while min_degree < max_degree:
    degree = min_degree + (max_degree - min_degree)//2
    print("Try %d" % degree)
    p = polynomial_through(roots=roots, degree=degree, xs=xs)
    if p is None:
        min_degree = degree + 1
    else:
        max_degree = degree

p = polynomial_through(roots=roots, degree=min_degree, xs=xs, verbose=True)
#p = (xs[0] - 1)**2 + (xs[1] - 1)**2 - 2
fp = sympy.lambdify(xs, p, modules='numpy')
print("Polynomium:")
print(sympy.Eq(p, 0))
print(sympy.factor(p, deep=True))
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
