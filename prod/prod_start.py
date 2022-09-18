import numpy as np
import glob

frames = [320,340,360,380,400,420,440,460,480,500]
for i,frame in enumerate(frames):
    for filename in glob.glob('../traj/CB1_assym_eq_1.nc'):
        f = open('../cpp/CB1_assym_pr_1_frame_'+str(i),'w')
        f.write('parm ../CB1_assym.prmtop \n')
        f.write('trajin ' + filename + '\n')
        f.write('autoimage origin \n')
        f.write('center origin \n')
        f.write('trajout ' + '../rst/CB1_assym_pr_1_frame_'+str(i) + '.rst7 onlyframes ' + str(frame) + '\n')
        f.write('run \n')
        f.write('quit \n')
        f.close()
