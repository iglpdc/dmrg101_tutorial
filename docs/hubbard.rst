DMRG algorithm for the Hubbard model
====================================

The goal of this exercise is to implement the finite version of the DMRG
algorithm for the one-dimensional Hubbard model. The Hamiltonian is:

.. math::
    H=-t\sum_{i, \sigma}\left(c^{\dagger}_{i, \sigma}c_{i, \sigma} + h.c.\right)+
    U\sum_{i}n_{i, \uparrow}n_{i, \downarrow}
 
where :math:`c_{i, \sigma}` is the destruction operator for an electron at
site :math:`i` and spin :math:`\sigma`, and
:math:`n_{i, \sigma}=c^{\dagger}_{i, \sigma}c_{i, \sigma}`.

Exercise
--------

Decide by yourself what you want to calculate with this code.

Solution
--------

The implementation is straighforward now. You have to write a model class
which with the functions to set the Hamiltonian, block Hamiltonians, and
operators to update. A possible implementation is:

.. literalinclude:: ../solutions/hubbard_helpers.py
    :pyobject: HubbardModel

You also have to write a single site class for the
electronic site, providing operators defined in the Hilbert space of the
site:

.. literalinclude:: ../solutions/hubbard_helpers.py
    :pyobject: ElectronicSite

The implementation of the DMRG algorithm itself has only minor changes
with respect what you have done before:

.. literalinclude:: ../solutions/tfim.py
    :pyobject: main

See :download:`a full implementation of the above code
<../solutions/tfim.py>`. To learn how to run this code you
can use:
::
    
    $ ./solutions/tfim.py --help
    Implements the full DMRG algorithm for the S=1/2 TFIM.

    Calculates the ground state energy and wavefunction for the
    Ising model in a transverse field for a chain of spin one-half. The
    calculation of the ground state is done using the full DMRG algorithm,
    i.e. first the infinite algorithm, and then doing sweeps for
    convergence with the finite algorithm.

    Usage:
      tfim.py (-m=<states> -n=<sites> -s=<sweeps> -H=<field>) [--dir=DIR -o=FILE]
      tfim.py -h | --help

    Options:
      -h --help         Shows this screen.
      -n <sites>        Number of sites of the chain.
      -m <states>       Number of states kept.
      -s <sweeps>       Number of sweeps in the finite algorithm.
      -H <field>        Magnetic field in units of coupling between spins.
      -o --output=FILE  Ouput file [default: tfim.dat]
      --dir=DIR         Ouput directory [default: ./]

