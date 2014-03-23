import math

import matplotlib as mpl
import matplotlib.ticker

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

class LatexFormatter(mpl.ticker.FormatStrFormatter):
    def __init__(self):
        super(LatexFormatter, self).__init__('$%g$')

import matplotlib.pyplot as plt
import numpy as np

roots = (-1, 1, 4, 6)
xrange = (-2, 8)
yrange = (-40, 40)

values = (20, 0, -10, 0)

def plot_fun(ax, f, x, *args):
    ax.plot(x, f(x), *args)

plt.axes([0.125,0.2,0.95-0.125,0.95-0.2])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.xaxis.set_major_formatter(LatexFormatter())
ax.yaxis.set_major_formatter(LatexFormatter())
ax.axis(xrange + yrange)
ax.set_xticks(range(xrange[0], xrange[1]+1, 1))

plot_fun(ax, lambda x: np.product([(x - root) for root in roots], axis=0),
        np.linspace(xrange[0], xrange[1], 500), 'b')
ax.plot(roots, np.zeros_like(roots), 'bx')

ax.plot(roots, values, 'r+')

ax.legend((r'$p \in \mathrm{Pol}_4(\mathbb{R})$', r'$\mathrm{rest}_S(p)$', r'$v \in \mathbb{R}^S$'))

width_pt = 341
width_in = 341 / 72.27
golden_mean = (math.sqrt(5)-1.0)/2.0
fig.set_size_inches(width_in, width_in * golden_mean)
fig.savefig('rest.pgf')
