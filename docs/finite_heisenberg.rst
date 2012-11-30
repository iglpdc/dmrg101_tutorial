Finite DMRG algorithm for the Heisenberg chain
==============================================

The goal of this exercise is to implement the finite version of the DMRG
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
model for a chain of spins :math:`S=1/2` using the finite version of the
DMRG algorithm. You should be able to pass the number of sites of the
chain, the number of states kept during the DMRG truncation, and the
number of sweeps during the finite algorithm as parameters. The output
should be the energy per site, entanglement entropy and, truncation error
at each step of the algorithm. You could:

- make a plot of the energy per site versus system size to see how
  converges. The relative error in energy is given by the truncation
  error. In a given sweep, where is the energy best approximate to the
  actual value?
- calculate the energy per site for a given system sizes and different
  number of states kept, make a plot, and extrapolate the values of the
  energy to zero truncation error with a linear fit. Compare this with the
  exact value for a finite size given by the Bethe Anstatz.
- find the scaling of the entanglement entropy for the Heisenberg model.
  The entanglement entropy for open boundary conditions scales as
  :math:`S(x)=\frac{c}{3}sin\left(\frac{\pi x}{2L}\right)`, where
  :math:`c, x, L` are the central charge of the model; the size of the
  part of the chain you keep in the reduced density matrix, and the size
  of the whole chain, respectively.

Solution
--------

This is the first implementation of a full DMRG algorithm in the tutorial.
The full DMRG algorithm involves doing first the infinite version of the
DMRG algorithm, which we cover in the last exercise, and then a number of
sweeps keeping the size of the system fixed. The latter are called finite
algorithm sweeps. During the finite sweeps, one block is growing exactly
as during the infinite version of the algorithm, while the other has to
shrink to keep the system size constant. As the truncation in DMRG can
take you only for larger block sizes, you have to use an old version of
the block that is shrinking. This is done saving each of the block when
they grow, so you simply pull the right one from a list of old blocks.

The first change we will make is to grow the system asymmetrically. This
means that during the infinite version of the algorithm we keep on block
(the right one) one-site long. You do this simply but not growing the
right block:

.. literalinclude:: ../solutions/heisenberg.py
    :pyobject: infinite_dmrg_step

Next thing is to write the function to implement the DMRG step during the
finite algorithm. The only difference with the infinite version is that
now in addition to grow one of the blocks, you set the other one to be the
one with the proper size in a previous sweep.

.. literalinclude:: ../solutions/heisenberg.py
    :pyobject: finite_dmrg_step

During the finite sweeps you want to increase the number of states that
you keep. A good way to do that is increasing them linearly at each
half-sweep from the number of states at the end of the infinite algorithm
to the number you want to keep at the end.

The rest is just doing a loop for the sweeps, each finite sweep comprising
a "half-sweep" to the left and a "half-sweep" to the right, and being sure
that during the finite algorithm the size of the system is constant. The
implementation (with some extras for saving the results and the rest)
looks like this:

.. literalinclude:: ../solutions/heisenberg.py
    :pyobject: main

See :download:`a full implementation of the above code
<../solutions/heisenberg.py>`. To learn how to run this code you
can use:
::
    
    $ ./tutorial/solutions/heisenberg.py --help

------------------------------------------------------------------------------
