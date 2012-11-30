Infinite DMRG algorithm for the Heisenberg chain
================================================

The goal of this exercise is to implement the infinite version of the DMRG
algorithm for the antiferromagnetic Heisenberg chain of spins one-half. 

The Hamiltonian that we will build is the antiferromagnetic Heisenberg
model:

.. math::
    H=\sum_{i}\vec{S}_{i}\cdot\vec{S}_{i+1}=
    \sum_{i}\left[S^{z}_{i}S^{z}_{i+1}+
    \frac{1}{2}\left(S^{\dagger}_{i}S^{-}_{i+1}+
    S^{-}_{i}S^{\dagger}_{i+1}\right)\right]

 
where :math:`\vec{S}_{i}=\vec{\sigma}_{i}/2` are the spin operators,
:math:`\sigma_{i}` are the Pauli matrices, and
:math:`S^{\pm}_{i}=S^{x}_{i}\pm i S^{y}_{i}`.

Exercise
--------

Calculate the ground state energy of the antiferromagnetic Heisenberg
model for a chain of spins :math:`S=1/2` using the infinite version of
the DMRG algorithm. You should be able to pass the number of sites of the
chain and the number of states kept during the DMRG truncation as
parameters. The output should be the energy per site, entanglement entropy
and, truncation error at each step of the algorithm. You could:

- make a plot of the energy per site versus system size to see how
  converges. The relative error in energy is given by the truncation
  error.
- calculate the energy per site for a given system sizes and different
  number of states kept, make a plot, and extrapolate the values of the
  energy to zero truncation error with a linear fit. Compare this with the
  exact value for a finite size given by the Bethe Anstatz.
- extract the central charge for the Heisenberg model. The entanglement
  entropy for open boundary conditions scales as 
  :math:`S(x)=\frac{c}{3}sin\left(\frac{\pi x}{2L}\right)`, where :math:`c,
  x, L` are the central charge of the model; the size of the part of the
  chain you keep in the reduced density matrix, and the size of the whole
  chain, respectively. In the infinite version of DMRG you are restricted to
  :math:`x=L/2`, but it enough to make a guess.

Solution
--------

You can use the `System` class to create the four spin chain, and set the
system Hamiltonian to the antiferromagnetic Heisenberg model, as we did in
the previous exercise. Then you calculate ground state as before. From
there you need to implement the DMRG transformation. We did most of the
required steps (calculating, diagonalizing, and truncating the reduced
density matrix) in a previous exercise. 

The new thing is updating the operators using the transformation matrix.
You need to specify two things before doing that. First which are going to
be the operators to be transformed. These are the ones that you need in
the next step of the DMRG algorithm. To build the Heisenberg Hamiltonian
using the block-site-site-block construction, you need your block to have
a truncated version of the current single site spin operators
:math:`S^{z}, S^{\dagger}, S^{-}`. Additionally you need to update the
block hamiltonian, which belongs to the current block and contains the
part of the Hamiltonian involving only degrees of freedom *within* the
block. This function sets the operators to be updated:

.. literalinclude:: ../solutions/infinite_heisenberg.py
    :pyobject: set_operators_to_update_to_AF_Heisenberg

The block hamiltonian is built by including the previous block hamiltonian, 
if any, and the terms coming from the interactions between the block and
the single site that is being included in the new block:

.. literalinclude:: ../solutions/infinite_heisenberg.py
    :pyobject: set_block_hamiltonian_to_AF_Heisenberg

Now we can put together the function to grow a block by one site:

.. literalinclude:: ../solutions/infinite_heisenberg.py
    :pyobject: grow_block_by_one_site

A step of the infinite DMRG algorithm consists in growing at the same time
both the left and the right blocks by one site. 

.. literalinclude:: ../solutions/infinite_heisenberg.py
    :pyobject: infinite_dmrg_step

The other thing remaining is to make a small change in the function to set
up the AF Heisenberg Hamiltonian that we used in the previous exercise to
include the block hamiltonians:

.. literalinclude:: ../solutions/infinite_heisenberg.py
    :pyobject: set_hamiltonian_to_AF_Heisenberg

Now you just make a loop that repeats the DMRG step until you reach the
desired system size. Adding some code to pass the arguments from the
command line and to save the results to a file, you get something like
this:

.. literalinclude:: ../solutions/infinite_heisenberg.py
    :pyobject: main

See :download:`a full implementation of the above code
<../solutions/infinite_heisenberg.py>`. To learn how to run this code you
can use:
::
    
    $ ./tutorial/solutions/infinite_heisenberg.py --help
    Implements the infinite version of the DMRG algorithm for the S=1/2 AF
    Heisenberg.
    
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

------------------------------------------------------------------------------
