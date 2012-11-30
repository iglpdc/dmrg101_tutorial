#!/usr/bin/env python
"""Implements the full DMRG algorithm for the S=1/2 AF Heisenberg.

Calculates the ground state energy and wavefunction for the
antiferromagnetic Heisenberg model for a chain of spin one-half. The
calculation of the ground state is done using the full DMRG algorithm,
i.e. first the infinite algorithm, and then doing sweeps for convergence
with the finite algorithm.

Usage:
  heisenberg.py (-m=<states> -n=<sites> -s=<sweeps>) [--dir=DIR -o=FILE]
  heisenberg.py -h | --help

Options:
  -h --help         Shows this screen.
  -n <sites>        Number of sites of the chain.
  -m <states>       Number of states kept.
  -s <sweeps>       Number of sweeps in the finite algorithm.
  -o --output=FILE  Ouput file [default: heisenberg.dat]
  --dir=DIR         Ouput directory [default: ./]

"""
from dmrg101.core.calculate_states_to_keep import calculate_states_to_keep
from dmrg101.core.entropies import calculate_entropy, calculate_renyi
from dmrg101.core.reduced_DM import diagonalize, truncate
from dmrg101.core.sites import SpinOneHalfSite
from dmrg101.core.system import System
from dmrg101.core.truncation_error import calculate_truncation_error
import dmrg101.utils.models.heisenberg_model as model
from docopt import docopt
import os

def grow_block_by_one_site(growing_block, ground_state_wf, system, 
		           number_of_states_kept):
    """Grows one side of the system by one site.

    Calculates the truncation matrix by calculating the reduced density
    matrix for `ground_state_wf` by tracing out the degrees of freedom of
    the shrinking side. Then updates the operators you need in the next
    steps, effectively growing the size of the block by one site. 	
    
    Parameters
    ----------
    growing_block : a string.
        The block which is growing. It must be 'left' or 'right'.
    ground_state_wf : a Wavefunction.
        The ground state wavefunction of your system.
    system : a System object.
        The system you want to do the calculation on. This function
	assumes that you have set the Hamiltonian to something.
    number_of_states_kept : an int.
        The number of states you want to keep in each block after the
	truncation. If the `number_of_states_kept` is smaller than the
	dimension of the current Hilbert space block, all states are kept.
 
    Returns
    -------
    entropy : a double.
        The Von Neumann entropy for the cut that splits the chain into two
	equal halves.
    truncation_error : a double.
        The truncation error, i.e. the sum of the discarded eigenvalues of
	the reduced density matrix.
    """
    system.set_growing_side(growing_block)
    rho = ground_state_wf.build_reduced_density_matrix(system.shrinking_side)
    evals, evecs = diagonalize(rho)
    truncated_evals, truncation_matrix = truncate(evals, evecs,
		                                  number_of_states_kept)
    entropy = calculate_entropy(truncated_evals)
    truncation_error = calculate_truncation_error(truncated_evals)
    model.set_block_hamiltonian(system)
    model.set_operators_to_update(system)
    system.update_all_operators(truncation_matrix)
    return entropy, truncation_error

def infinite_dmrg_step(system, number_of_states_kept):
    """Performs one step of the (asymmetric) infinite DMRG algorithm.

    Calculates the ground state of a system with a given size, then
    performs the DMRG transformation on the operators of *one* block,
    therefore increasing by one site the number of sites encoded in the
    Hilbert space of this block, and reset the block in the system to be
    the new, enlarged, truncated ones. The other block is kept one-site
    long.

    Parameters
    ----------
    system : a System object.
        The system you want to do the calculation on. This function
	assumes that you have set the Hamiltonian to something.
    number_of_states_kept : an int.
        The number of states you want to keep in each block after the
	truncation. If the `number_of_states_kept` is smaller than the
	dimension of the current Hilbert space block, all states are kept.
 
    Returns
    -------
    energy : a double.
        The energy for the `current_size`.
    entropy : a double.
        The Von Neumann entropy for the cut that splits the chain into two
	equal halves.
    truncation_error : a double.
        The truncation error, i.e. the sum of the discarded eigenvalues of
	the reduced density matrix.

    Notes
    -----
    This asymmetric version of the algorithm when you just grow one of the
    block while keeping the other one-site long, is obviously less precise
    than the symmetric version when you grow both sides. However as we are
    going to sweep next using the finite algorithm we don't care much
    about precision at this stage.
    """
    model.set_hamiltonian(system)
    ground_state_energy, ground_state_wf = system.calculate_ground_state()
    entropy, truncation_error = grow_block_by_one_site('left', ground_state_wf, 
		                                       system, number_of_states_kept)
    return ground_state_energy, entropy, truncation_error

def finite_dmrg_step(growing_block, system, left_block_size, number_of_states_kept):
    """Performs one step of the finite DMRG algorithm.

    Calculates the ground state of a system with a given size, then
    performs the DMRG transformation on the operators of *one* block,
    therefore increasing by one site the number of sites encoded in the
    Hilbert space of this block, and reset the block in the system to be
    the new, enlarged, truncated ones. The other block is read out from
    the previous sweep.

    Parameters
    ----------
    growing_block : a string.
        The block which is growing. It must be 'left' or 'right'.
    system : a System object.
        The system you want to do the calculation on. This function
	assumes that you have set the Hamiltonian to something.
    left_block_size : an int.
        The number of sites in the left block in the *current* step, not
	including the single site.     
    number_of_states_kept : an int.
        The number of states you want to keep in each block after the
	truncation. If the `number_of_states_kept` is smaller than the
	dimension of the current Hilbert space block, all states are kept.
 
    Returns
    -------
    energy : a double.
        The energy at this step.
    entropy : a double.
        The Von Neumann entropy for the cut at this step.
    truncation_error : a double.
        The truncation error, i.e. the sum of the discarded eigenvalues of
	the reduced density matrix.

    Raises
    ------
    DMRGException
        if `growing_side` is not 'left' or 'right'.

    Notes
    -----
    This asymmetric version of the algorithm when you just grow one of the
    block while keeping the other one-site long, is obviously less precise
    than the symmetric version when you grow both sides. However as we are
    going to sweep next using the finite algorithm we don't care much
    about precision at this stage.
    """
    model.set_hamiltonian(system)
    ground_state_energy, ground_state_wf = system.calculate_ground_state()
    if growing_block not in ('left', 'right'):
	raise DMRGException('Growing side must be left or right.')

    entropy, truncation_error = grow_block_by_one_site(growing_block, 
		                                       ground_state_wf,
		                                       system, number_of_states_kept)
    system.set_block_to_old_version(left_block_size)
    return ground_state_energy, entropy, truncation_error

def main(args):
    # 
    # create a system object with spin one-half sites and blocks.
    #
    spin_one_half_site = SpinOneHalfSite()
    system = System(spin_one_half_site)
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
	    infinite_dmrg_step(system, 
		               number_of_states_infinite_algorithm) )
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
                finite_dmrg_step('right', system, left_block_size, states) )
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
                finite_dmrg_step('left', system, left_block_size, states) )
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
