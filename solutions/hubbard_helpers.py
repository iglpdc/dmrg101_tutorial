""" A file with some stuff to use in the docs.

This file is **not** meant to be used in the actual code. It just copies a
couple of functions from the dmrg101.core package to display them in the
tutorial docs.
"""

class ElectronicSite(Site):
    """A site for electronic models
    
    You use this site for models where the single sites are electron
    sites. The Hilbert space is ordered such as:

    - the first state, labelled 0,  is the empty site,
    - the second, labelled 1, is spin down, 
    - the third, labelled 2, is spin up, and 
    - the fourth, labelled 3, is double occupancy.
    
    Notes
    -----
    Postcond: The site has already built-in the spin operators for: 

    - c_up : destroys an spin up electron,
    - c_up_dag, creates an spin up electron,
    - c_down, destroys an spin down electron,
    - c_down_dag, creates an spin down electron,
    - s_z, component z of spin,
    - s_p, raises the component z of spin,
    - s_m, lowers the component z of spin,
    - n_up, number of electrons with spin up,
    - n_down, number of electrons with spin down,
    - n, number of electrons, i.e. n_up+n_down, and
    - u, number of double occupancies, i.e. n_up*n_down.

    """
    def __init__(self):
        super(ElectronicSite, self).__init__(4)
	# add the operators
        self.add_operator("c_up")
        self.add_operator("c_up_dag")
        self.add_operator("c_down")
        self.add_operator("c_down_dag")
        self.add_operator("s_z")
        self.add_operator("s_p")
        self.add_operator("s_m")
        self.add_operator("n_up")
        self.add_operator("n_down")
        self.add_operator("n")
        self.add_operator("u")
	# for clarity
        c_up = self.operators["c_up"]
        c_up_dag = self.operators["c_up_dag"]
        c_down = self.operators["c_down"]
        c_down_dag = self.operators["c_down_dag"]
        s_z = self.operators["s_z"]
        s_p = self.operators["s_p"]
        s_m = self.operators["s_m"]
        n_up = self.operators["n_up"]
        n_down = self.operators["n_down"]
        n = self.operators["n"]
        u = self.operators["u"]
	# set the matrix elements different from zero to the right values
	# TODO: missing s_p, s_m
        c_up[0,2] = 1.0
        c_up[1,3] = 1.0
        c_up_dag[2,0] = 1.0
        c_up_dag[3,1] = 1.0
        c_down[0,1] = 1.0
        c_down[2,3] = 1.0
        c_down_dag[1,0] = 1.0
        c_down_dag[3,2] = 1.0
        s_z[1,1] = -1.0
        s_z[2,2] = 1.0
        n_up[2,2] = 1.0
        n_up[3,3] = 1.0
        n_down[1,1] = 1.0
        n_down[3,3] = 1.0
        n[1,1] = 1.0
        n[2,2] = 1.0
        n[3,3] = 2.0
        u[3,3] = 1.0

class HubbardModel(object):
    """Implements a few convenience functions for Hubbard model.
    
    Does exactly that.
    """
    def __init__(self):
        super(HubbardModel, self).__init__()
		
    def set_hamiltonian(self, system):
        """Sets a system Hamiltonian to the Hubbard Hamiltonian.
    
        Does exactly this. If the system hamiltonian has some other terms on
        it, there are not touched. So be sure to use this function only in
        newly created `System` objects.
    
        Parameters
        ----------
        system : a System.
            The System you want to set the Hamiltonian for.
        """
        system.clear_hamiltonian()
        if 'bh' in system.left_block.operators.keys():
            system.add_to_hamiltonian(left_block_op='bh')
        if 'bh' in system.right_block.operators.keys():
            system.add_to_hamiltonian(right_block_op='bh')
        system.add_to_hamiltonian('c_up', 'c_up_dag', 'id', 'id', -1.)
        system.add_to_hamiltonian('c_up_dag', 'c_up', 'id', 'id', -1.)
        system.add_to_hamiltonian('c_down', 'c_down_dag', 'id', 'id', -1.)
        system.add_to_hamiltonian('c_down_dag', 'c_down', 'id', 'id', -1.)
        system.add_to_hamiltonian('id', 'c_up', 'c_up_dag', 'id', -1.)
        system.add_to_hamiltonian('id', 'c_up_dag', 'c_up', 'id', -1.)
        system.add_to_hamiltonian('id', 'c_down', 'c_down_dag', 'id', -1.)
        system.add_to_hamiltonian('id', 'c_down_dag', 'c_down', 'id', -1.)
        system.add_to_hamiltonian('id', 'id', 'c_up', 'c_up_dag', -1.)
        system.add_to_hamiltonian('id', 'id', 'c_up_dag', 'c_up', -1.)
        system.add_to_hamiltonian('id', 'id', 'c_down', 'c_down_dag', -1.)
        system.add_to_hamiltonian('id', 'id', 'c_down_dag', 'c_down', -1.)
        system.add_to_hamiltonian('u', 'id', 'id', 'id', self.U)
        system.add_to_hamiltonian('id', 'u', 'id', 'id', self.U)
        system.add_to_hamiltonian('id', 'id', 'u', 'id', self.U)
        system.add_to_hamiltonian('id', 'id', 'id', 'u', self.U)
    
    def set_block_hamiltonian(self, system):
        """Sets the block Hamiltonian to the Hubbard model block Hamiltonian.
    
        Parameters
        ----------
        system : a System.
            The System you want to set the Hamiltonian for.
        """
        # If you have a block hamiltonian in your block, add it
        if 'bh' in system.growing_block.operators.keys():
            system.add_to_block_hamiltonian('bh', 'id')
        system.add_to_block_hamiltonian('c_up', 'c_up_dag', -1.)
        system.add_to_block_hamiltonian('c_up_dag', 'c_up', -1.)
        system.add_to_block_hamiltonian('c_down', 'c_down_dag', -1.)
        system.add_to_block_hamiltonian('c_down_dag', 'c_down', -1.)
        system.add_to_block_hamiltonian('id', 'u', self.U)
        system.add_to_block_hamiltonian('u', 'id', self.U)
    
    def set_operators_to_update(self, system):
        """Sets the operators to update to the ones for the Hubbard model.
    
        Parameters
        ----------
        system : a System.
            The System you want to set the Hamiltonian for.
        """
        # If you have a block hamiltonian in your block, update it
        if 'bh' in system.growing_block.operators.keys():
            system.add_to_operators_to_update('bh', block_op='bh')
        system.add_to_operators_to_update('c_up', site_op='c_up')
        system.add_to_operators_to_update('c_up_dag', site_op='c_up_dag')
        system.add_to_operators_to_downdate('c_down', site_op='c_down')
        system.add_to_operators_to_downdate('c_down_dag', site_op='c_down_dag')
        system.add_to_operators_to_update('u', site_op='u')
