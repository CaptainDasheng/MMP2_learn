#==============================================================================
#==============================================================================

# contributor
__contributor__ = 'Dongwen Yang'

#revised date
__date__ = Apr. 4, 2017 



class SetHSE(object):

     
    """
    class SetHSE :: Prepare the input for the HSE calculations, mainly
          modfiy the kpoints.
    """
    
    def __init__(self, vasp, lkpoints=10, args, **kwargs):


        self.functional = vasp


    # initalizing the kpoints % 
    def init_kmesh(self):

        """
        :: initalizing the kmesh within whole BZ,
        
        return str
        """

    # interpolation of kpoints % 
    def linspace_kpoints(self, hpoints, lkpoints=10, *args, **kwargs):

        """
        :: interpolation along the direction of high symmetry path.

              hpoints: a couple of kpoints to be considered, such as, 
                 tuple({'\Gamma':[0,0,0]}, {'L':[0.5,0.5,0.5]})

              lkpoints: the number of the kpoints along hpoints, 
                 default is 10.

              interkpoints: the final coordination of interpolated k.

        return lkpoints, str(interkpoints)

        """

        pass

    # organize the kpoints % 
    def __orgnized_kpoints__(self):

        """
        modfiy the kpoints the functional.

        return kpoints 
        """
        
        kpoints 

        self.functional.kpoints = kpoints 
