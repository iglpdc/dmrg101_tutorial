#!/usr/bin/env python2.7
"""Implements the infinite version of the DMRG algorithm for the S=1/2 AF Heisenberg.

Calculates the ground state energy and wavefunction for the
antiferromagnetic Heisenberg model for a chain of spin one-half. The
calculation of the ground state is done using the infinite version of the
DMRG algorithm.

Usage:
  infinite_heisenberg.py (-m=<states> -n=<sites>) [--dir=DIR -o=FILE]
  infinite_heisenberg.py -h | --help

Options:
  -h --help         Shows this screen.
  -n <sites>        Number of sites of the chain.
  -m <states>       Number of states kept.
  -o --output=FILE  Ouput file [default: infinite_heisenberg.dat]
  --dir=DIR         Ouput directory [default: ./]

"""
from docopt import docopt
import sys, os
sys.path.insert(0, os.path.abspath('../../dmrg101'))
from dmrg101.core.entropies import calculate_entropy, calculate_renyi
from dmrg101.core.reduced_DM import diagonalize, truncate
from dmrg101.core.sites import SpinOneHalfSite
from dmrg101.core.system import System

def set_hamiltonian_to_AF_Heisenberg(system):
    """Sets a system Hamiltonian to the AF Heisenberg Hamiltonian.

    Does exactly this. If the system hamiltonian has some other terms on
    it, there are not touched. So be sure to use this function only in
    newly created `System` objects.

    Parameters
    ----------
    system : a System.
        The System you want to set the Hamiltonain for.
    """
    system.clear_hamiltonian()
    if 'bh' in system.left_block.operators.keys():
        system.add_to_hamiltonian(left_block_op='bh')
    if 'bh' in system.right_block.operators.keys():
        system.add_to_hamiltonian(right_block_op='bh')
    system.add_to_hamiltonian('id', 'id', 's_z', 's_z')
    system.add_to_hamiltonian('id', 'id', 's_p', 's_m', .5)
    system.add_to_hamiltonian('id', 'id', 's_m', 's_p', .5)
    system.add_to_hamiltonian('id', 's_z', 's_z', 'id')
    system.add_to_hamiltonian('id', 's_p', 's_m', 'id', .5)
    system.add_to_hamiltonian('id', 's_m', 's_p', 'id', .5)
    system.add_to_hamiltonian('s_z', 's_z', 'id', 'id')
    system.add_to_hamiltonian('s_p', 's_m', 'id', 'id', .5)
    system.add_to_hamiltonian('s_m', 's_p', 'id', 'id', .5)

def set_block_hamiltonian_to_AF_Heisenberg(system):
    """Sets the block Hamiltonian to be what you need for AF Heisenberg.

    Parameters
    ----------
    system : a System.
        The System you want to set the Hamiltonian for.
    """
    # If you have a block hamiltonian in your block, add it
    if 'bh' in system.growing_block.operators.keys():
        system.add_to_block_hamiltonian('bh', 'id')
    system.add_to_block_hamiltonian('s_z', 's_z')
    system.add_to_block_hamiltonian('s_p', 's_m', .5)
    system.add_to_block_hamiltonian('s_m', 's_p', .5)

def set_operators_to_update_to_AF_Heisenberg(system):
    """Sets the operators to update to be what you need to AF Heisenberg.

    Parameters
    ----------
    system : a System.
        The System you want to set the Hamiltonian for.
    """
    # If you have a block hamiltonian in your block, update it
    if 'bh' in system.growing_block.operators.keys():
        system.add_to_operators_to_update('bh', block_op='bh')
    system.add_to_operators_to_update('s_z', site_op='s_z')
    system.add_to_operators_to_update('s_p', site_op='s_p')
    system.add_to_operators_to_update('s_m', site_op='s_m')

def infinite_dmrg_step(system, current_size, number_of_states_kept):
    """Performs one step of the infinite DMRG algorithm.

    Calculates the ground state of a system with a given size, then
    performs the DMRG transformation on the operators of *both* blocks,
    therefore increasing by one site the number of sites encoded in the
    Hilbert space of each blocks, and reset the blocks in the system to be
    the new, enlarged, truncated ones.

    Parameters
    ----------
    system : a System object.
        The system you want to do the calculation on. This function
	assumes that you have set the Hamiltonian to something.
    current_size : an int.
        The number of sites in the chain. It must be even and it is a
	constant in this function, only used to calculate the energy per
	site.
    number_of_states_kept : an int.
        The number of states you want to keep in each block after the
	truncation. If the `number_of_states_kept` is smaller than the
	dimension of the current Hilbert space block, all states are kept.
 
    Returns
    -------
    current_size : an int.
        The same parameter that you passed.
    energy_per_site : a double.
        The energy per site for the `current_size`.
    entropy_left : a double.
        The Von Neumann entropy for the cut that splits the chain into two
	equal halves.

    Notes
    -----
    Normally you don't update both blocks. If the chain is symmetric, you
    just can use the operators for the one of the sides to mirror the
    operators in the other side, saving the half of the CPU time. In
    practical DMRG calculations one uses the finite algorithm to
    improve the result of the infinite algorithm, and one of the blocks
    is kept one site long, and therefore not updated.
    """
    set_hamiltonian_to_AF_Heisenberg(system)
    ground_state_energy, ground_state_wf = system.calculate_ground_state()
    # do the left block
    rho = ground_state_wf.build_reduced_density_matrix('left')
    evals, evecs = diagonalize(rho)
    truncated_evals, truncation_matrix = truncate(evals, evecs,
		                                  number_of_states_kept)
    entropy_left = calculate_entropy(truncated_evals)
    system.set_growing_side('left')
    set_block_hamiltonian_to_AF_Heisenberg(system)
    set_operators_to_update_to_AF_Heisenberg(system)
    system.update_all_operators(truncation_matrix)
    # do the right block
    rho = ground_state_wf.build_reduced_density_matrix('right')
    evals, evecs = diagonalize(rho)
    truncated_evals, truncation_matrix = truncate(evals, evecs,
		                                  number_of_states_kept)
    entropy_right = calculate_entropy(truncated_evals)
    system.set_growing_side('right')
    set_block_hamiltonian_to_AF_Heisenberg(system)
    set_operators_to_update_to_AF_Heisenberg(system)
    system.update_all_operators(truncation_matrix)
    # as a check you could  print `entropy_right`, it should be equal to
    # `entropy_left`
    # print entropy_right
    return current_size, ground_state_energy/current_size, entropy_left

def main(args):
    # 
    # create a system object with spin one-half sites and blocks.
    #
    spin_one_half_site = SpinOneHalfSite()
    system = System(spin_one_half_site)
    number_of_sites = int(args['-n'])
    number_of_states_kept= int(args['-m'])
    #
    # infinite DMRG algorithm
    #
    number_of_sites = 2*(number_of_sites/2) # make it even
    for current_size in range(4, number_of_sites+1, 2):
	current_size, energy_per_site, entropy = (
	    infinite_dmrg_step(system, current_size, number_of_states_kept) )
	print current_size, current_size*energy_per_site, entropy
	current_size += 2

if __name__ == '__main__':
    args = docopt(__doc__, version = 0.1)
    main(args)
    output_file = os.path.join(os.path.abspath(args['--dir']), args['--output'])
    print 'Results stored in ' + output_file
