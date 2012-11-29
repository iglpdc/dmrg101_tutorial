#!/usr/bin/env python2.7
""" Calculates the entanglement entropy of a two qbit system

Calculates the von Neumann entanglement entropu of a system of two
spin one-half spins restricted to the subspace of total spin equal to
zero. 

Usage:
  two_qbit_system.py [--dir=DIR -o=FILE]
  two_qbit_system.py -h | --help

Options:
  -h --help         Shows this screen.
  -o --output=FILE  Ouput file [default: two_qbit_entropies.dat]
  --dir=DIR         Ouput directory [default: ./]

"""
from docopt import docopt
import sys, os
sys.path.insert(0, os.path.abspath('../../dmrg101'))
from math import cos, sin, pi
from dmrg101.core.entropies import calculate_entropy
from dmrg101.core.reduced_DM import diagonalize
from dmrg101.core.wavefunction import Wavefunction

def create_two_qbit_system_in_singlet(psi):
    """ Returns the wf of the system as a function of `psi`.

    The (normalized) wavefunction of the two-qbit system can be
    parametrized as a function an angle `psi`. 

    Parameters
    ----------
    psi : a double 
        Parametrizes the wavefunction.
    
    Returns
    -------
    result : a Wavefunction
        The wavefunction of the two-qbit system for the given `psi`.
    """
    result = Wavefunction(2, 2)

    # set the different components.
    result.as_matrix[0, 0] = 0.
    result.as_matrix[0, 1] = cos(psi)
    result.as_matrix[1, 0] = sin(psi)
    result.as_matrix[1, 1] = 0.
    return result

def trace_out_left_qbit_and_calculate_entropy(wf):
    """Calculates the entropy after tracing out the left qbit.

    To calculate the entanglement entropy you need to first build the
    reduced density matrix tracing out the degrees of freedom of one of
    the two qbits (it does not matter which, we pick up left here.)

    Parameters
    ----------
    wf : a Wavefunction
        The wavefunction you build up the reduced density matrix with.

    Returns
    -------
    result : a double
        The value of the von Neumann entanglement entropy after tracing
	out the left qbit.
    """
    reduced_DM_for_right_qbit = wf.build_reduced_density_matrix('left')
    evals, evecs = diagonalize(reduced_DM_for_right_qbit)
    result = calculate_entropy(evals)
    return result

def main(args):
    """Calculates the entanglement entropy for a system of two qbits in a
    singlet state.
    """
    # 
    # get a bunch of values (number_of_psi) for psi
    #
    number_of_psi = 1000
    step = 2*pi/number_of_psi
    psi_values = [x*step for x in range(number_of_psi)] 
    #
    # python function map applies a function to a sequence
    #
    wfs = map(create_two_qbit_system_in_singlet, psi_values)
    entropies = map(trace_out_left_qbit_and_calculate_entropy, wfs)
    # 
    # find to which value of psi corresponds the max entropy
    #
    zipped = zip(psi_values, entropies)
    max_value = max(zipped, key=lambda item: (item[1]))
    # 
    # print the results
    #
    print "The maximum value for entropy is %8.6f." %max_value[1]
    print "The wavefunction with max entropy is: "
    print create_two_qbit_system_in_singlet(max_value[0]).as_matrix
    #
    # save for plotting
    #
    filename = args['--output']
    f = open(filename, 'w')
    f.write('\n'.join('%s %s' % x for x in zipped))
    f.close()

if __name__ == '__main__':
    args = docopt(__doc__, version = 0.1)
    main(args)
    output_file = os.path.join(os.path.abspath(args['--dir']), args['--output'])
    print "The whole list of psi vs entropies is saved in",
    print output_file+'.'
