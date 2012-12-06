DMRG algorithm for the Ising model in a transverse field
========================================================

The goal of this exercise is to implement the finite version of the DMRG
algorithm for the Ising model in a transverse field (TFIM) for a chain of
spins one-half and study the model close to the quantum phase transition.

The Hamiltonian is:

.. math::
    H=\sum_{i}\left(-JS^{z}_{i}S^{z}_{i+1}-hS^{x}_{i}\right)

 
where :math:`\vec{S}_{i}=\vec{\sigma}_{i}/2` are the spin operators, and
:math:`\sigma_{i}` are the Pauli matrices.

This model has a quantum phase transition at :math:`h/J=1.0`. At the
critical point the exact value for the ground state energy for a finite
system with open boundary conditions is given by:

.. math::
    E/J=1-cosec\left(\frac{\pi}{2\left(2L+1\right)}\right)

The high field phase :math:`h/J>1.0` is a paramagnet with an order
parameter :math:`\langle S^{x}\rangle\neq 0`. The low field phase
:math:`h/J<1.0` is a ferromagnet with an order parameter :math:`\langle
S^{z}\rangle\neq 0`. At the critical point, both order parameters go to
zero.

The exact values of the energies at the critical point are:

+----+---------------+
| L  | E/L (J)       |
+====+===============+
| 16 | -1.2510242438 |
+----+---------------+
| 20 | -1.255389856  |
+----+---------------+
| 32 | -1.2620097863 |
+----+---------------+
| 40 | -1.264235845  |
+----+---------------+
| 64 | -1.267593439  |
+----+---------------+

Exercise
--------

Study the quantum phase transition in the TFIM by measuring the order
parameters at each side of the transition and comparing the energy at the
critical point with the exact result. The critical point is 

- calculate the energy per site for a given system size and different
  number of states kept, make a plot, and extrapolate the values of the
  energy to zero truncation error with a linear fit. Compare this with the
  exact value for a finite size given by the exact solution.
- measure the value of the order parameters for a given system size and
  a few values of :math:`h/J` around the critical point. Use a reasonable
  system size and a number of states, so you actually are able to get to
  sample the phase diagram around the critical point. Get the order
  parameter by summing up the values of :math:`S^{x}` or :math:`S^{z}`
  for all sites of the chain. Plot both order parameters versus :math:`h/J`.

Solution
--------

The implementation goes is pretty similar to the one for the Heisenberg
model of the last exercise. The main change is of course the change of the
Hamiltonian, block Hamiltonian, and the operators you need to update after
each DMRG step. 

To simplify things and switching between different models, the `System`
class has a few convience functions to do grow the blocks, perform the
infinite and finite DMRG steps, and set the Hamiltonian. These are the
same functions we have seen before, but now written as methods of the
`System` class. 

When using these methods you have to write a `Model` class that contains
the details of the TFIM model:

.. literalinclude:: ../static/tfim_helpers.py
    :pyobject: TranverseFieldIsingModel

The best thing is to look at the final implementation and compare to the
previous one for the Heisenberg model:

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

