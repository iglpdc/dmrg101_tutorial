#!/usr/bin/env python
""" Plots data from a text file using matplotlib.

The file must contain two columns of the same lenght. The first column
will be the data in the x-axis, the second column the data in the y-axis.

Usage:
  two_qbit_plot.py <file>
  two_qbit_plot.py -h | --help

Options:
  -h --help       Shows this screen.

"""
from docopt import docopt
import matplotlib.pyplot as plt
import os

args = docopt(__doc__, version = 0.1)
data_file = os.path.abspath(args['<file>'])
fin = open(data_file, 'r')
x = []
y = []
for lines in fin:
    line = lines.split()
    x.append(float(line[0]))
    y.append(float(line[1]))
plt.plot(x, y)
plt.show()
