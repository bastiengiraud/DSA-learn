# dynrun.py - run this script from PSSe GUI
import os, sys

# SET VARs
study   = 'ieee118_flat'
savkey  = 'ieee118'
dyr     = 'ieee118'
dyradd  = ''
cnvfile = '%s_cnv'%savkey
snpfile = '%s.snp'%savkey
outfile = '%s.out'%study
logfile = '%s.log'%study
# -------------------------------------------------------------------------
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
_s = psspy.getdefaultchar()

# SET LOG FILE [uncomment next two lines to save progress output to file]
psspy.progress_output(2,logfile,[0,0])
# -------------------------------------------------------------------------
# 1: set solution_parameters
# -------------------------------------------------------------------------
psspy.solution_parameters_3([_i,100,_i],
                            [_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
# -------------------------------------------------------------------------
# 2: convert case and create snp file
# -------------------------------------------------------------------------
#** convert loads and create converted case
if not os.path.isfile(cnvfile):
    psspy.case(savkey)
    psspy.runrspnsfile(r"""Conl.idv""")
    psspy.cong()
    psspy.ordr(0)
    psspy.fact()
    psspy.tysl(0)
    psspy.save(cnvfile)
else:
    psspy.case(cnvfile)
#** Read DYRE records -  system dynamics + solar PV dynamics
if not os.path.isfile(snpfile):
    psspy.dyre_new([1,1,1,1],
                    '%s.dyr'%dyr,
                    r"""conec.for""",
                    r"""conet.for""","")
    if dyradd:
        psspy.dyre_add([_i,_i,_i,_i],
                        '%s.dyr'%dyradd,
                        r"""conecb.for""",
                        r"""conetb.for""")
    
    #** Save snapshot for dynamics
    psspy.snap([-1,-1,-1,-1,-1],snpfile)
else:
    psspy.rstr(snpfile)
# -------------------------------------------------------------------------
# 3: run 
#    No user-define models were called, therefore no need to compile
#    add channels
#    select a large gen to set the relative angle option
#    update solution parametes as applicable
# -------------------------------------------------------------------------
#** Set channels
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,2,0])
psspy.chsb(0,1,[-1,-1,-1,1,3,0])
psspy.chsb(0,1,[-1,-1,-1,1,4,0])

#** Run dynamics simulation
psspy.set_relang(1, 89, '1')
psspy.dynamics_solution_params([99,_i,_i,_i,_i,_i,_i,_i],
                               [0.8,_f, 0.001, 0.004,_f,_f,_f,_f],'')
psspy.set_chnfil_type(0)    #1 for OUTX format, 0 for (old) OUT format 
psspy.strt(1,outfile)

#apply a disturbance at t= 1.0, normal clear after 5 cycles
#  [to test under no disturbance, comment next five lines]
#psspy.run(0,1.0,99,9,0)
#psspy.dist_bus_fault(77,1,1.0,[0.0,-0.2E+10])
#psspy.run(0,1.08,99,3,0)            #LVRT test = 9 cycles under fault
#psspy.dist_clear_fault(1)
#psspy.dist_branch_trip(77,82,r"""1""")

#continue with run
psspy.run(0,10.0,501,9,0)
#psspy.snap([-1,-1,-1,-1,-1],snpfile)
psspy.progress_output(1)
