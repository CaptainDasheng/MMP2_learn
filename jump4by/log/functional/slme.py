#!/usr/bin/env python 

class CalcSLME(object):

    """
    class CalcSLME aiming to obtain the maximium theoretical solar
          cell efficiency based on the dielectric matrix.
            
    params:
            
            stdout:: running path of vasp optical calculation, default is
                     current path;

            thickness:: the thickness of films, default is 3.0 miu meter;

            dielect:: file include the dielectric matrix, default is
                     OUTCAR; alternative file is vasprun.xml;

            shifts:: scissor to modify the band gap to the experiemental
                     values or accuracy values obtained via more accuracy 
                     functionals, such HSE/GW, default is zero;
                    
    """

        def __init__(self, stdout=None, thickness=3.0, output='OUTCAR', shifts=0.0):
            
            

            if stdout is None:

                path = os.getcwd()

            if output is 'OUTCAR':

                self.dielectric = get_epsilon_outcar()

            elif dielect is 'vasprun.xml':

                self.dielectric = get_epsilon_vasprun()
            
            else:
		# other program under devolped %
                raise IOError('invalid the input file 
                              contained dielectric matrix.')

            self.calc_slme()


        def get_epsilon_outcar(self):

            pass

        def get_epsilon_vasprun(self):

            pass

        def slme(self):
            
            am15 = open("AM15-eta_SQ.dat",'r').readlines() # solar spectrum 
            ev = 1.60217648740E-19
            h = 6.626068E-34
            c = 299792458
            Pin = 1000/ev


            slme_val  = []
            indirect  = self.indirect
            direct    = self.direct 
            shifts    = self.shifts
            dielect_matrix  = self.dielectric_matrix

            thickness = np.linspace(0, self.thickness, 1000)
            
            for t in thickness:

                s = self.cal_slme(dielect_matrix, t*e-4, indirect, direct, shifts)


        def calc_slme(self, dielect_matrix, ):




                     
    
