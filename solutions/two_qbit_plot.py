#!/usr/bin/env python
""" Plots data from a text file using matplotlib.

The file must contain two columns of the same lenght. The first column
will be the data in the x-axis, the second column the data in the y-axis.

Usage:
  two_qbit_system.py [--dir=DIR -o=FILE]
  two_qbit_system.py -h | --help

Options:
  -h --help       Shows this screen.
  -f --file=FILE  Ouput file [default: two_qbit_entropies.dat]
  --dir=DIR        Ouput directory [default: ./]

"""
from docopt import docopt
import matplotlib.pyplot as plt

args = docopt(__doc__, version = 0.1)
data_file = os.path.join(os.path.abspath(args['--dir']), args['--file'])
fin = open(data_file, 'r')
x = []
y = []
for lines in fin:
    line = lines.split()
    x.append(float(line[0]))
    y.append(float(line[1]))
plt.plot(x, y)
plt.show()
