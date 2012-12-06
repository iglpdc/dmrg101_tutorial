#!/usr/bin/env python
""" Plots data from a text file using matplotlib.

The file must contain two columns of the same lenght. The first column
will be the data in the x-axis, the second column the data in the y-axis.

Usage:
  plot_from_file.py <file> [-x=<x_col> -y=<y_col>]
  plot_from_file.py -h | --help

Options:
  -h --help   Shows this screen.
  -x <x_col>  Which column has the data for x-axis [default: 0]
  -y <y_col>  Which column has the data for y-axis [default: 1]

"""
from docopt import docopt
import matplotlib.pyplot as plt
import os

args = docopt(__doc__, version = 0.1)
data_file = os.path.abspath(args['<file>'])
x_col = int(args['-x'])
y_col = int(args['-y'])
fin = open(data_file, 'r')
x = []
y = []
for lines in fin:
    line = lines.split()
    x.append(float(line[x_col]))
    y.append(float(line[y_col]))
plt.plot(x, y)
plt.show()
