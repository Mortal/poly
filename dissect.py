import itertools
import numpy as np
import sympy
import random

def itermonomials_degree(xs, degree):
    if not degree:
        yield sympy.S.One
        return
    x, *xs = xs
    if not xs:
        yield x ** degree
        return
    yield x ** degree
    for d in range(1, degree):
        for m in itermonomials_degree(xs, d):
            yield x ** (degree - d) * m
    yield from itermonomials_degree(xs, degree)

def itermonomials(xs):
    """Yield all monomials over `xs` in grlex order."""
    for degree in itertools.count():
        yield from itermonomials_degree(xs, degree)

def coefficient_tuples(Xs, monomials):
    for m in monomials:
        d = m.as_powers_dict()
        yield tuple(int(d[X]) for X in Xs)

def coefficient_names(tuples):
    for m in tuples:
        yield 'c%s' % (''.join('_%s' % exponent for exponent in m))

def coefficients(Xs, monomials):
    return sympy.symbols(list(coefficient_names(coefficient_tuples(Xs, monomials))))

def polynomial_through(roots, Xs):
    k = len(roots)
    ms = list(itertools.islice(itermonomials(Xs), 1+k))
    system = sympy.Matrix([[m.subs(zip(Xs, x)) for m in ms] + [0] for x in roots])
    symbols = coefficients(Xs, ms)
    if k == 1:
        return Xs[0] - roots[0][0]
    sol = sympy.solvers.minsolve_linear_system(system, *symbols, quick=False)
    pol = sum(sol[s] * m for s, m in zip(symbols, ms) if s in sol)
    #print("Solution:")
    #print(pol)
    #for k, v in sol.items():
    #    print("  %s = %s" % (k, v))
    return pol

def dissecting_polynomial_trial(As, Xs, threshold):
    #print(' '.join(str(A.shape) for A in As))
    ai = [random.choice(A) for A in As]
    pol = polynomial_through(ai, Xs)
    f = sympy.lambdify(Xs, pol, modules='numpy')
    well_dissecting = 0
    good = []
    bad = []
    for A in As:
        t = threshold*len(A)
        values = f(*A.T)
        above = (values > 0)
        below = (values < 0)
        #print(above)
        #print(below)
        if above.sum() <= t and below.sum() <= t:
            good.append((A[above, :], A[below, :]))
        else:
            bad.append(A)
    return good, bad, pol, f

def unit_vectors(d):
    return sympy.symbols(list('x%d' % (1+i) for i in range(d)))

def well_dissecting_polynomial(As, threshold=7/8):
    As = tuple(np.array(A) for A in As)

    # k: number of point sets
    k = len(As)

    try:
        ni, di = zip(*[Ai.shape for Ai in As])
    except ValueError:
        raise ValueError("input should be a list of matrices with points as rows")
    #print("ni = %s" % list(ni))

    # n: total number of points
    n = sum(ni)
    #print("n = %d points in total" % n)

    # d: dimension of containing space
    ds = frozenset(di)
    try:
        d, = ds
    except ValueError:
        raise ValueError("different dimensions in input: %s" % sorted(ds))

    #print("Dimension is %d" % d)

    Xs = unit_vectors(d)
    well_dissecting = 0
    while well_dissecting < (1+len(As))//2:
        good, bad, pol, f = dissecting_polynomial_trial(As, Xs, threshold)
        well_dissecting = len(good)
    return good, bad, pol

def partitioning_polynomial(P, r):
    P = np.array(list(P))
    n, d = P.shape
    print("r-partitioning n points in dimension d")
    print("    r = %s" % r)
    print("    n = %d" % n)
    print("    d = %d" % d)
    Xs = unit_vectors(d)
    polynomials = []
    sets = [P]
    product = sympy.S.One
    j = 0
    def print_row(**kwargs):
        print("{j}\t{s}\t{p}\t{l}".format(**kwargs))

    print_row(j='j', s='s', p='|Pj|', l='|Lj(s)|')

    print_row(j=0, s='', p=len(sets), l='')

    while max(len(each) for each in sets) > n/r:
        j += 1
        threshold = int((7/8)**j * n)
        large = tuple(filter(lambda each: len(each) > threshold, sets))
        sets = list(filter(lambda each: len(each) <= threshold, sets))
        if not large:
            continue

        #print("Phase %d threshold: %d" % (j, threshold))
        #print("Large set sizes: %s" % [len(each) for each in large])
        #print("Small set sizes: %s" % [len(each) for each in sets])

        remaining = tuple(large)
        #f = sympy.S.One
        gs = []
        s = 0
        while remaining:
            print_row(j='' if s else j, s=s, p=len(sets), l=len(remaining))
            s += 1
            good, bad, g = well_dissecting_polynomial(remaining)
            #print("In phase (%d, %d), well dissected %d, leaving %d" %
            #        (j, s, len(good), len(bad)))
            # Flatten `good`
            for above, below in good:
                sets.append(above)
                sets.append(below)
            remaining = tuple(bad)
            gs.append(g)
            #f *= g
        print_row(j='', s=s, p=len(sets), l=len(remaining))
        polynomials.append(gs)
        #product *= f
        #print("Product is %s" % product)
    return polynomials, sets

def main():
    n = 2000
    d = 3
    A = (n * np.random.random((n, d))).astype(np.int)
    r = 50
    #ni = 10
    #As = tuple(zip(*[((0, a), (a, 0), (ni, a), (a, ni)) for a in range(ni)]))
    #print(well_dissecting_polynomial(As))
    #p = partitioning_polynomial(itertools.chain.from_iterable(As), 4)
    p, sets = partitioning_polynomial(A, r)
    print(',\n'.join('[%s]' % ',\n '.join(map(str, factors)) for factors in p))
    print("Expecting sets of size at most %s" % (len(A)/r))
    print("Set sizes: %s" % sorted(len(each) for each in sets))
    print("Exceptional points: %s" % (len(A) - sum(len(each) for each in sets)))
    #sols = sympy.solve(p, dict=True)
    #for sol in sols:
    #    print(sol)
    #    dim = d - len(sol)
    #    print("That's a %d-surface!" % dim)
    #    is_flat = True
    #    for dependent, expr in sol.items():
    #        if expr.has(*sol.keys()):
    #            raise Exception("Solution has interdependent variables")

    #        terms, gens = expr.as_terms()
    #        for term, data in terms:
    #            coeff, monom, ncpart = data
    #            if ncpart:
    #                raise Exception("Term has non-commutative part")
    #            if monom and max(monom) > 1:
    #                is_flat = False
    #    if is_flat:
    #        print("It's flat!")


if __name__ == '__main__':
    main()
