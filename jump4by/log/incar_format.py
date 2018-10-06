incarf=\
"""
#=======================jump2 set incar=======================================
#  system name % 
   SYSTEM     = {name}        
          
#=============================================================================
# start parameter % 
 {tprecise}   PREC          = {precise}   # normal or accurate (medium, high low)
 {tistart}    ISTART        = {istart}    # job   : 0-new  1-cont  2-samecut
 {ticharg}    ICHARG        = {icharg}    # charge: 1-file 2-atom 10-const
 {tispin}     ISPIN         = {ispin}     # spin polarized calculation?
 {tlsorbit}   LSORBIT       = {lsorbit}   # spin-orbit coupling
 {tlasph}     LASPH         = {lasph}     # aspherical Exc in radial PAW
 {tmetagga}   METAGGA       = {metagga}   # non-selfconsistent MetaGGA calc.
 {tlnoncollinear} LNONCOLLINEAR = {lnoncollinear}     # non collinear calculations

#=============================================================================
# electronic relaxation % 
  {tencut}   ENCUT  = {encut}           # energy cutoff  
  {tnelm}    NELM   = {nelm}            # threshold step for electronic 
  {tediff}   EDIFF  = {ediff}           # stopping-criterion for ELM
  {tlreal}   LREAL  = {lreal}           # real-space projection

# optimize electronic relaxation % 
  {tiaglo} IALGO = {iaglo}    
  {tamix}  AMIX  = {amix}
  {tamin}  AMIN  = {amin}
 
# electrons and spin % 
  {tnelect} NELECT = {nelect}          # total number of electrons

#=============================================================================
# ionic relaxation % 
  {tediffg}   EDIFFG  = {ediffg}          # stopping-criterion for IOM
  {tnsw}      NSW     = {nsw}          # number of steps for IOM
  {tnblock}   NBLOCK  = {nblock}          # inner block;
  {tkblock}   KBLOCK  = {kblock}          # outer block 
  {tibrion}   IBRION  = {ibrion}          # ionic relax: 0-MD 1-quasi-New 2-CG
  {tisif}     ISIF    = {isif}          # stress and relaxation
  {tisym}     ISYM    = {isym}          # 0-nonsym 1-usesym 2-fastsym
             
  {tpotim}    POTIM   = {potim}          # time-step for ionic-motion
  {ttemin}    TEIN    = {temin}          # initial temperature
  {ttebeg}    TEBEG   = {tebeg}          # temperature during run
  {tteend}    TEEND   = {teend}  
  {tsmass}    SMASS   = {smass}          # Nose mass-parameter (am)
  {tpstress}  PSTRESS = {pstress}          # pullay stress

#=============================================================================
# DOS related values % 
  {temin}        EMIN    = {emin}           # energy-range for DOS
  {temax}        EMAX    = {emax}           
  {tismear}      ISMEAR  = {ismear}            # broadening -4-tet -1-fermi 0-gaus
  {tsigma}       SIGMA   = {sigma}           # 

#=============================================================================
# write flags % 
  {tlwave}  LWAVE   = {lwave}	             # write WAVECAR
  {tlcharg} LCHARG  = {lcharg}	             # write CHGCAR
  {tlvtot}  LVTOT   = {lvtot}	             # write LOCPOT, total local potential
  {tlvhar}  LVHAR   = {lvhar}	             # write LOCPOT, Hartree potential only
  {tlelf}   LELF    = {lelf}	             # write ELF
  {tlorbit} LORBIT  = {lorbit}	             # 0 simple, 1 ext, 2 COOP (PROOUT)

#=============================================================================
# exchange correlation treatment % 
  {tgga}      GGA     = {gga}	             # GGA type
  {tlhfcalc}  LHFCALC = {lhfcalc}	             # Hartree Fock is set to
  {taexx}     AEXX    = {aexx}               # exact exchange contribution
                                       
#=============================================================================
# linear response parameters % 
  {tlepsilon}  LEPSILON = {lepsilon}              # dielectric tensor
  {tlrpa}      LRPA     = {lrpa}              # RPA
  {tcshift}    CSHIFT   = {cshift}         # complex shift 
                                       
#=============================================================================
# dimension of arrays %      # 
  {tnbands} NBANDS   = {nbands}	     # number of bands    
  {tnedos}  NEDOS    = {nedos}	     # number of dos      

  {tng} NGX  = {ngx};    NGY  = {ngy};    NGZ  = {ngz}   # dimension x,y,z 
  {tngf} NGXF = {ngxf};   NGYF = {ngyf};   NGZF = {ngzf}  # dimension x,y,z 
#=============================================================================
# parallelisation %
  {tnpa} NPAR   = {npar}
  {tlplane} LPLANE = {lplane}
"""   
