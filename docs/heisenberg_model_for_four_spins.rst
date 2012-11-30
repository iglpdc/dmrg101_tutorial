Heisenberg model for four spins
===============================

The goal of this exercise is to build a Hamiltonian for the four spin
system, and calculate its ground state. The difference with the previous
exercise is that we will use the same setup as for the DMRG algorithm.
This means that we will build a chain as a block-site-site-block system,
where the left-most spin will be the left block, the two central spins
will be the single sites, and the right-most spin will be the right block.

Again, the Hamiltonian that we will build is the antiferromagnetic Heisenberg
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

Calculate the ground state (energy and wavefunction) of the
antiferromagnetic Heisenberg model for a system of four spins one-half.

Solution
--------

You can use the `System` class to create the four spin chain, and set the
system Hamiltonian to the antiferromagnetic Heisenberg model.  To do that
you add the terms to the system Hamiltonian, and to add a term you just
have to specify which operator acts in each of the parts: left block, left
site, right site, and right block, using the name of the operator. You can
specify a parameter if you need to with a fifth optional argument. This
function does that:

.. literalinclude:: ../solutions/heisenberg_for_four_spins.py
    :pyobject: set_hamiltonian_to_AF_Heisenberg

The remaining is calling again to the Lanczos solver to get the ground
state energy and wavefunction. You can use a convenience function in the
`System` class. Putting everything together you end up with something like
this:

.. literalinclude:: ../solutions/heisenberg_for_four_spins.py
    :pyobject: main

See :download:`a full implementation of the above code
<../solutions/heisenberg_for_four_spins.py>`. If you run that code you should
get:
::
    
    (dmrg101) $ python tutorial/solutions/heisenberg_for_four_spins.py
    The ground state energy is -1.616025.
    The ground state wavefunction is: 
    [[ ...

------------------------------------------------------------------------------
