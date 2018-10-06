#===========================================
#
#
#===========================================



class CalcBandStructure(object):

    """
    Aiming the calculate the band structure via scf and nonscf 
    calculations. 
    
    three types of band structures can be calculated:
        a). band structure with line_mode;
        b). band structure with line_mode by using high precise methods,
            e.g. hybrid functional and green wavefunction;
        c). band structure with in 3D space;

    methods:
        func:: the object of vasp;
        high_symmetry_path:: get the high symmetry path;
        scf_band:: self-consistent calculation;
        nonscf_band:: nonself-consistent calculation of band structures;
        calc_bandflow:: the workflow to calculate band structures;

    """

    def __init__(self, func, *args, **kwargs):
        
        if isinstance(func, Vasp):
           self.func = func 

        self.calc_bandflow()

    def calc_bandflow(self):
        
        func = self.func 
        
        # get primary cell % 
        if func.to_primary is True:
            func.structure = func.structure.get_primary_cell()

        # self-consistent calculation % 
        nonscf = func.calculator('band/scf')
        nonscf.icharg = 11
        nonscf.istart = 1
        nonscf.nband  = nonscf.get_nband()*2.0
        for path in nonscf.high_symmetry_points:
            nonscf.caclulator(path)


