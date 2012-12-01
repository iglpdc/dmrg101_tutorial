#!/usr/bin/env python
"""Implements the full DMRG algorithm for the Hubbard model.

Calculates the ground state energy and wavefunction for the
Ising model in a transverse field for a chain of spin one-half. The
calculation of the ground state is done using the full DMRG algorithm,
i.e. first the infinite algorithm, and then doing sweeps for convergence
with the finite algorithm.

Usage:
  tfim.py (-m=<states> -n=<sites> -s=<sweeps> -H=<U_over_t>) [--dir=DIR -o=FILE]
  tfim.py -h | --help

Options:
  -h --help         Shows this screen.
  -n <sites>        Number of sites of the chain.
  -m <states>       Number of states kept.
  -s <sweeps>       Number of sweeps in the finite algorithm.
  -U <U_over_t>     Electronic interaction in units of hopping.
  -o --output=FILE  Ouput file [default: hubbard.dat]
  --dir=DIR         Ouput directory [default: ./]

"""
from dmrg101.core.calculate_states_to_keep import calculate_states_to_keep
from dmrg101.core.sites import ElectronicSite 
from dmrg101.core.system import System
from dmrg101.utils.models.hubbard_model import HubbardModel
from docopt import docopt
import os

def main(args):
    # 
    # create a system object with electron sites and blocks, and set
    # its model to be the Hubbard model.
    #
    electronic_site = ElectronicSite()
    system = System(electronic_site)
    system.model = HubbardModel()
    #
    # read command-line arguments and initialize some stuff
    #
    number_of_sites = int(args['-n'])
    number_of_states_kept= int(args['-m'])
    number_of_sweeps= int(args['-s'])
    number_of_states_infinite_algorithm = 10
    if number_of_states_kept < number_of_states_infinite_algorithm:
	number_of_states_kept = number_of_states_infinite_algorithm
    sizes = []
    energies = []
    entropies = []
    truncation_errors = []
    system.number_of_sites = number_of_sites
    #
    # infinite DMRG algorithm
    #
    max_left_block_size = number_of_sites - 3
    for left_block_size in range(1, max_left_block_size):
	energy, entropy, truncation_error = ( 
	    system.infinite_dmrg_step(number_of_states_infinite_algorithm) )
	current_size = left_block_size + 3
	sizes.append(left_block_size)
	energies.append(energy/current_size)
	entropies.append(entropy)
	truncation_errors.append(truncation_error)
    #
    # finite DMRG algorithm
    #
    states_to_keep = calculate_states_to_keep(number_of_states_infinite_algorithm, 
		                              number_of_states_kept,
		                              number_of_sweeps)
    half_sweep = 0
    while half_sweep < len(states_to_keep):
	# sweep to the left
        for left_block_size in range(max_left_block_size, 1, -1):
	    states = states_to_keep[half_sweep]
            energy, entropy, truncation_error = ( 
                system.finite_dmrg_step('right', left_block_size, states) )
            sizes.append(left_block_size)
            energies.append(energy/number_of_sites)
            entropies.append(entropy)
            truncation_errors.append(truncation_error)
	half_sweep += 1
	# sweep to the right
	# if this is the last sweep, stop at the middle
	if half_sweep == 2 * number_of_sweeps - 1:
	    max_left_block_size = number_of_sites / 2
        for left_block_size in range(1, max_left_block_size):
            energy, entropy, truncation_error = ( 
                system.finite_dmrg_step('left', left_block_size, states) )
            sizes.append(left_block_size)
            energies.append(energy/number_of_sites)
            entropies.append(entropy)
            truncation_errors.append(truncation_error)
	half_sweep += 1
    # 
    # save results
    #
    output_file = os.path.join(os.path.abspath(args['--dir']), args['--output'])
    f = open(output_file, 'w')
    zipped = zip (sizes, energies, entropies, truncation_errors)
    f.write('\n'.join('%s %s %s %s' % x for x in zipped))
    f.close()
    print 'Results stored in ' + output_file

if __name__ == '__main__':
    args = docopt(__doc__, version = 0.1)
    main(args)
