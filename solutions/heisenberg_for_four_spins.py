#!/usr/bin/env python2.7
""" Calculates the ground state for the AF Heisenberg for two spins 1/2.

Calculates the ground state energy and wavefunction for the
antiferromagnetic Heisenberg model for a chain of four spin one-half. The
calculation of the ground state is done using the Lanczos algorithm.

Usage:

  heisenberg_for_four_spins.py -h | --help

Options:
  -h --help         Shows this screen.

"""
from docopt import docopt
import sys, os
sys.path.insert(0, os.path.abspath('../../dmrg101'))
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
    system.add_to_hamiltonian('id', 'id', 's_z', 's_z')
    system.add_to_hamiltonian('id', 'id', 's_p', 's_m', .5)
    system.add_to_hamiltonian('id', 'id', 's_m', 's_p', .5)
    system.add_to_hamiltonian('id', 's_z', 's_z', 'id')
    system.add_to_hamiltonian('id', 's_p', 's_m', 'id', .5)
    system.add_to_hamiltonian('id', 's_m', 's_p', 'id', .5)
    system.add_to_hamiltonian('s_z', 's_z', 'id', 'id')
    system.add_to_hamiltonian('s_p', 's_m', 'id', 'id', .5)
    system.add_to_hamiltonian('s_m', 's_p', 'id', 'id', .5)

def main():
    # 
    # create a system object with spin one-half sites and blocks.
    #
    spin_one_half_site = SpinOneHalfSite()
    system = System(spin_one_half_site)
    #
    # build the Hamiltonian, and solve it using Lanczos.
    #
    set_hamiltonian_to_AF_Heisenberg(system)
    ground_state_energy, ground_state_wf = system.calculate_ground_state()
    #
    # print results
    #
    print "The ground state energy is %8.6f." %ground_state_energy
    print "The ground state wavefunction is :"
    print ground_state_wf.as_matrix

if __name__ == '__main__':
    args = docopt(__doc__, version = 0.1)
    main()
