Playing with a two-qbit system
==============================

The goal of this exercise is to build the wavefunction of a pair of spins
one-half (a.k.a. a pair of qbits), and calculate its entanglement
entropy. This is a pretty trivial exercise that you could do without much
hassle in a piece of paper (see below.) The purpose is just get you
familiar with how things are done in the code before moving to bigger
adventures.

Exercise 
-------- 

Calculate the entanglement entropy when you trace out one of the spins in
a general state in the subspace of a two spin one-half system where the
two particles have opposite spin. For which of the states in this
subspace is are the two spins maximally entangled (i.e. the Von Neumann
entanglement entropy is maximal)?

Hint
----

The first thing we need to do is to write the wavefunction of the two spin
system. In dmrg101 the wavefunctions are represented as matrices instead
of vectors, which may be more familiar to you. 

The reason for that is that as in DMRG we always have to split the
physical systems (say a chain of spins, or the two spins of the problem)
in left and right subsystems, the notation with matrices is more suited.
In the dmrg101 code the rows of the matrix representing a wavefunction
correspond to states of the left subsystem, and the columns correspond to
states of the right subsystem.

For example to represent the two spin one-half system, with one spin as
left subsystem and the other as right subsystem, we need a 2x2 matrix.
Matrix elements in the first row (column) will correspond to states with
the left (right) spin down.  Matrix elements in the second row (column)
will correspond to states with the left (right) spin up. The choice of
whether the first or second row corresponds to spin down or up is
arbitrary, but once you made the choice you have to be consistent.

If we restrict ourselves to the subspace where the particles have
opposite spin, the most general wavefunction for the two qbit systems
is simply [#]_:

.. math::
    |\psi\rangle = \cos \phi |\downarrow\uparrow\rangle 
    + \sin \phi |\uparrow\downarrow\rangle = 
    \begin{pmatrix} 0 & \cos \phi \\ \sin\phi & 0 \end{pmatrix}

Solution
--------

The plan is the following. First we are going to write a function to
calculate the wavefunction for the two-qbit system as a function of 
an angle `psi`:

.. literalinclude:: ../solutions/two_qbit_system.py
    :pyobject: create_two_qbit_system_in_singlet

Now we are going to get the reduced density matrix tracing out the
left qbit and calculate the corresponding entanglement entropy:

.. literalinclude:: ../solutions/two_qbit_system.py
    :pyobject: trace_out_left_qbit_and_calculate_entropy

Now it just a matter to generate a bunch of different values for `psi`,
calculate the corresponding wavefunction with the first function above,
and pass the wavefunction to the second funciton above to get the value
for the entropy. The following code makes this:

.. literalinclude:: ../solutions/two_qbit_system.py
    :pyobject: main

See :download:`a full implementation of the above code
<../solutions/two_qbit_system.py>`. If you run that code you should
get something like this: 
::
    
    $ ./solutions/two_qbit_system.py
    The maximum value for entropy is 0.693147.
    The wavefunction with max entropy is: 
    [[ 0.          0.70710678]
     [ 0.70710678  0.        ]]
    The whole list of psi vs entropies is saved in two_qbit_entropies.dat.

which are in fact the entropy (:math:`log(2)`) and the wavefunction of
the triplet state where the particles have opposite spins. In your own
code, you might alternatively have observed a state with the same entropy
but where one of the two components is negative; this is the singlet
state. See also the :download:`data for the entropies vs psi
<../solutions/two_qbit_entropies.dat>`. To see the options avaliable:
::

    $ ./solutions/two_qbit_system.py --help
    Calculates the entanglement entropy of a two qbit system
    
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

Conclusion
----------

It is important that you note that this is the general solution for a
system of two qbits, and that two-qbits cannot be more entangled than
when they are in an equal amplitude (but possibly opposite phase)
superposition state over configurations where the particles have opposite
spin. In system of many particles is splitted in two parts (think in a
larger chain of spins cut at some point in two), one can always represent
the relevant degrees of freedom at the cut as a set of qbits. Then it
follows from the result you just proved that the most *economical* way of
representing the entanglement across the cut is to map the degrees of
freedom of each side to a qbits and *maximally entangle* them across
the cut. Any other state to be formed with the qbits in one side and the
other, will either have less entanglement across the cut than the one
in the original degrees of freedom, or use more qbits at each side of
the cut. This is the basis of the mappings used in quantum information
methods like MPS or TNS, and you will see *maximally entangled
spins/qbits* a lot in the rest of the school.

--------------------------------------

.. [#] There could be an extra phase between the two components, but it
    cancels out later, so we don't bother to include it.

