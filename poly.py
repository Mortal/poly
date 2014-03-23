import sys
import matplotlib as mpl
if __name__ == '__main__' and '20points.pgf' in sys.argv:
    mpl.use('pgf')
    with open('matplotlib.tex') as f:
        mpl.rcParams.update({
            'font.family': 'serif',
            'axes.unicode_minus': False,
            'pgf.rcfonts': False,
            'pgf.preamble': f.read().splitlines(),

            'font.size': 10,
            'axes.labelsize': 10,
            'font.size': 10,
            'text.fontsize': 10,
            'legend.fontsize': 10,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
        })

import matplotlib.pyplot as plt
import sympy
import sympy.solvers
import numpy as np
import scipy
import scipy.ndimage
import math
import time
import itertools
import functools

def print_matrix(mat):
    if hasattr(mat, 'tolist'):
        mat = mat.tolist()
    rows = [list(map(str, row)) for row in mat]
    cols = tuple(zip(*rows))
    colwidth = max(max(map(len, row)) for row in rows)
    screen_width = 80
    indent = 2 * ' '
    space = 2 * ' '
    columns = (screen_width - len(indent) + len(space)) // (colwidth + len(space))
    for i in range(0, len(cols), columns):
        j = min(len(cols), i + columns)
        if i + 1 == j:
            print("Column %d:" % j)
        else:
            print("Columns %d - %d:" % (i+1, j))
        print("")
        for row in zip(*cols[i:j]):
            print('%s%s' % (indent, space.join(map(lambda s: s.rjust(colwidth), row))))
        print("")

def default_dimensions(n):
    """Default names (in order) given to axes in n dimensional space."""
    if n <= 3:
        names = ('x', 'y', 'z')
        return sympy.symbols(names[:n])
    else:
        return sympy.symbols(tuple('x%d' % n for n in range(1, n+1)))

def default_coefficients(n):
    """Default names (in order) given to entries of n-element coefficient vector."""
    if n <= 3:
        names = ('a', 'b', 'c')
        return sympy.symbols(names[:n])
    else:
        return sympy.symbols(tuple('a%d' % n for n in range(1, n+1)))

def monomials(xs, degree):
    """Yield all monomials over `xs` of degree at most `degree`."""
    return sorted(sympy.itermonomials(xs, degree),
            key=sympy.polys.orderings.monomial_key('grlex', tuple(reversed(xs))))

def coefficients(xs, degree):
    symbols = []
    e_str = '%d' if degree < 10 else '_%d'
    for p in monomials(xs, degree):
        m = sympy.Monomial(p, gens=xs)
        symbols.append('c%s' % ''.join(e_str % e for e in m.exponents))
    return tuple(sympy.symbols(symbols))

def remove_denominator(e):
    n, d = e.as_numer_denom()
    if d.is_constant():
        # Get rid of denominator
        return n
    else:
        print("Warning: Common denominator is")
        print(d)
        return e

def polynomial_through(roots, degree, xs, *, general=False, quick=True, verbose=False):
    ms = monomials(xs, degree)
    system = sympy.Matrix([[m.subs(zip(xs, x)) for m in ms] + [0] for x in roots])
    symbols = tuple(coefficients(xs, degree))
    print("%d roots, %d coefficients" % (len(roots), len(symbols)))

    if verbose:
        print("Linear system:")
        print_matrix(system)
        print(symbols)
    if general or not quick:
        sol = sympy.solvers.solve_linear_system(system, *symbols)
    else:
        sol = sympy.solvers.minsolve_linear_system(system, *symbols, quick=quick)
    if verbose:
        print("Solution:")
        for k, v in sol.items():
            print("  %s = %s" % (k, v))

    if () == tuple(filter(lambda n: n != 0, sol.values())):
        if verbose:
            print("No non-trivial solution")
        return None

    p = sum(sol.get(ci, ci) * m for ci, m in zip(symbols, ms))
    frees_orig = tuple(frozenset(p.free_symbols).difference(xs))
    frees = default_coefficients(len(frees_orig))
    p = p.subs(zip(frees_orig, frees))

    if verbose:
        print("Degrees of freedom: %d" % len(frees))
        print("General solution:")
        print(p)

    if frees and not general:
        # Choose concrete values for free variables
        # (Does not happen when minsolve_linear_system is used)
        concretes = []
        for free in frees:
            c = p.subs((v, 1 if v == free else 0) for v in frees)
            c = remove_denominator(c)
            concretes.append(c)
        degrees = map(lambda p: p.as_poly(*xs).total_degree(), concretes)
        ps, ds = zip(*sorted(zip(concretes, degrees), key=lambda a: a[1]))
        if verbose:
            print("Concretizations:")
            for pi, di in zip(ps, ds):
                print("  %s  (degree %d)" % (sympy.factor(pi), di))
        p = ps[0]

    return p

def min_degree_polynomial_through(roots, xs, *, final_verbose=True, **kwargs):
    max_degree = next(filter(lambda d: len(monomials(xs, d)) > len(roots), itertools.count(1)))
    min_degree = 1
    #min_degree = max_degree
    while min_degree < max_degree:
        print("Search for smallest feasible degree in [%d, %d]"
                % (min_degree, max_degree))
        # Don't do binary search here since running time grows fast with degree
        degree = min_degree + (max_degree - min_degree)//2
        degree = min_degree
        t1 = time.time()
        p = polynomial_through(roots=tuple(tuple(map(float, root)) for root in roots),
                degree=degree, xs=xs,
                **dict(tuple(kwargs.items()) + (('general', True),)))
        t2 = time.time()
        if p is None:
            print("Degree %d is infeasible (took %g)" % (degree, t2-t1))
            min_degree = degree + 1
        else:
            print("Degree %d is feasible (took %g)" % (degree, t2-t1))
            max_degree = degree

    print("Best degree is %d" % min_degree)

    return polynomial_through(roots=roots, degree=min_degree, xs=xs, quick=False,
            **dict(tuple(kwargs.items()) + (('verbose', final_verbose),)))

def check_roots(roots, fp):
    for xi in roots:
        print("Check root %s:" % (xi,))
        print(fp(*xi))

def plot_solution(ax, roots, fp, margin=1):
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

    x1 -= margin*w
    x2 += margin*w
    y1 -= margin*h
    y2 += margin*h

    xs = np.linspace(x1, x2, 500)
    ys = np.linspace(y1, y2, 500)
    xy = np.meshgrid(xs, ys)
    ax.contour(xy[0], xy[1], fp(xy[0], xy[1]), (0,))
    rootx, rooty = zip(*roots)
    ax.plot(rootx, rooty, 'o')

def make_plot(ax, dims, roots, **kwargs):
    xs = default_dimensions(dims)

    p = min_degree_polynomial_through(roots=roots, xs=xs, final_verbose=True)
    #p = (xs[0] - 1)**2 + (xs[1] - 1)**2 - 2

    print("Curve:")
    print(sympy.Eq(p, 0))

    fp = sympy.lambdify(xs, p, modules='numpy')

    check_roots(roots, fp)

    if dims == 2:
        plot_solution(ax, roots, fp, **kwargs)
    else:
        print("Not plotting since dimension is %d, not 2" % dims)

class LatexFormatter(mpl.ticker.FormatStrFormatter):
    def __init__(self):
        super(LatexFormatter, self).__init__('$%g$')

def main():
    if '20points.pgf' in sys.argv:
        dims = 2
        n = 20
        roots = tuple(map(tuple, 100*np.random.random((n, dims))))

        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.xaxis.set_major_formatter(LatexFormatter())
        ax.yaxis.set_major_formatter(LatexFormatter())

        make_plot(ax=ax, dims=dims, roots=roots, margin=0.1)

        width_pt = 341
        width_in = 341 / 72.27
        golden_mean = (math.sqrt(5)-1.0)/2.0
        fig.set_size_inches(width_in, width_in * golden_mean)
        fig.savefig('20points.pgf')
    elif 'time' in sys.argv:
        dims = 2
        xs = default_dimensions(dims)
        for degree in itertools.count(7):
            n = len(monomials(xs, degree)) - 1
            roots = tuple(map(tuple, 100*np.random.random((n, dims))))
            #degree = next(filter(lambda d: len(monomials(xs, d)) > len(roots), itertools.count(1)))
            t1 = time.time()
            p = polynomial_through(roots=roots, degree=degree, xs=xs)
            t2 = time.time()
            print("%d\t%d\t%g" % (degree, n, t2-t1))
            if not p:
                print("Infeasible")
    elif sys.argv[1] == 'image':
        dims = 2

        img_array = scipy.ndimage.imread(sys.argv[2], flatten=True)
        h, w = img_array.shape
        cs = scipy.cumsum(img_array.reshape(-1))
        cs /= cs[-1]
        n = 44
        vs = np.random.random(n)
        i = np.searchsorted(cs, vs)
        print(i)
        i, j = (i % w), (i // w)

        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.axis('equal')
        #ax.imshow(img_array)

        roots = list(zip(i, j))
        print(roots)
        make_plot(ax=ax, dims=dims, roots=roots)

        plt.show()
    else:

        dims = 2

        #roots = ((0, 0), (2, sympy.Rational(1,2)), (1, 2), (3, 3))

        #roots = tuple((2*i, 2+i) for i in range(7))

        degree = 12
        xs = default_dimensions(dims)
        n = len(monomials(xs, degree)) - 1
        roots = tuple(map(tuple, 100*np.random.random((n, dims))))
        #roots = [(x, sympy.exp(x - 1)) for x in [sympy.Rational(i, n-1) for i in range(n)]]

        #dims = 1
        #roots = [[x] for x in (-1, 1, 4, 6)]

        #roots = ((0, 1), (2, 0), (4, 1), (4, 4), (2, 5), (0, 4))

        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.axis('equal')

        make_plot(ax=ax, dims=dims, roots=roots)

        plt.show()

if __name__ == '__main__':
    main()
