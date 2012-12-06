#!/usr/bin/env python
""" Plots entanglement entropy data from a text file using matplotlib.

The file must contain two columns of the same lenght. The first column are
the sizes of the growing block, not including the single site. The second
column is not used.  The second column are the entanglement entropies
calculated for this cut.  The number of sites must be even and the lower
branch, i.e. when number of sites in the growing side -not including the
single site- is odd, is used. This is the ouput format for the scripts
that implement the DMRG algorithm in this tutorial.

The quantity plotted in the x-axis is :math:`ln(2x')`, where
:math:`x'=(L/\pi)\sin{\pi x /L}`. The slope of this plot is proportional
to the central charge of the CFT:

.. :math: S_{vN} \approx \frac{c}{6}ln(x')

The plot show the data as blue points and a green line with the CFT result
above, as a guide to the eye (not a real fit to data.)

See e.g. A. B. Kallin et al, Phys. Rev. Lett. 103, 117203 (2009).

Usage:
  plot_entropies.py <file> (-n=<sites>) 
  plot_entropies.py -h | --help

Options:
  -h --help       Shows this screen.
  -n <sites>        Number of sites of the chain.

"""
from docopt import docopt
import numpy as np
import math
import matplotlib.pyplot as plt
import os

args = docopt(__doc__, version = 0.1)
data_file = os.path.abspath(args['<file>'])
number_of_sites = int(args['-n'])
fin = open(data_file, 'r')
x = []
y = []
for lines in fin:
    line = lines.split()
    # the system (opposed to enviroment) has even sites:
    if int(line[0]) % 2 == 1:
	x_prime = math.sin(math.pi * (float(line[0]) + 1) / number_of_sites)
	x_prime *= 2 * number_of_sites / math.pi
        x.append(math.log(x_prime))
        y.append(float(line[2]))

# cft_result is the part of the CFT result proportional to the central
# charge (there are other smaller terms for finite number_of_sites),
# forced to pass by the last point of your data. So it's not a real fit,
# just a guide to the eye.
cft_result = np.array(x) / 6 
cft_result += y[-1] - x[-1] / 6
plt.plot(x, y, 'bo')
plt.plot(x, cft_result, 'g-')
plt.show()
