import glob
import numpy as np
import mdtraj as md
top = '../CB2_assym-stripped.prmtop'

for filename in glob.glob('./striped/CB2_assym_pr_*_frame_*'):
    t = md.load(filename,top=top)
    nframe= t.n_frames

    lig_1 = t.topology.select('(type C or type O or type N) and resid 319')
    lig_1_cor = np.mean(t.xyz[:,lig_1,:],axis=1)

    lig_2 = t.topology.select('(type C or type O or type N) and resid 320')
    lig_2_cor = np.mean(t.xyz[:,lig_2,:],axis=1)

    dist = np.empty([nframe,6])
    dist[:,0] = lig_1_cor[:,0]
    dist[:,1] = lig_1_cor[:,1]
    dist[:,2] = lig_1_cor[:,2]
    dist[:,3] = lig_2_cor[:,0]
    dist[:,4] = lig_2_cor[:,1]
    dist[:,5] = lig_2_cor[:,2]
    newfile = filename.split('-strip.dcd')[0].split('./striped/')[1]
    np.save('./numpy/'+newfile+'_lig_com',dist)
