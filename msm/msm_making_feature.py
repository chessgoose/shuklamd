import glob
import numpy as np
import mdtraj as md
from sklearn import preprocessing
import os
import sys
import pickle

if __name__=='__main__':
    top = '../CB2_assym-stripped.prmtop'
    totdist = []
    totfile = []
    totframes = []

    binding_pocket_residue = sorted([25,87,91,94,95,281,285,194,190,191,261,265,182,183,113,117,184])
    pairs = [[binding_pocket_residue[i],binding_pocket_residue[j]] for i in range(len(binding_pocket_residue)) for j in range(i+1,len(binding_pocket_residue)) if abs(binding_pocket_residue[i]-binding_pocket_residue[j]) > 4]

    for filename in glob.glob('./striped/CB2_assym_pr_*_frame_*'):
        #try:
        t = md.load(filename,top=top)
        nframe = t.n_frames
        dist  = np.empty([nframe,len(pairs)+6*len(binding_pocket_residue)])
        alpha_carbon = [[t.topology.select('name CA and resid ' + str(residue[0]-1))[0],t.topology.select('name CA and resid ' + str(residue[1]-1))[0]] for residue in pairs]
        dist[:,:len(pairs)] = md.compute_distances(t,alpha_carbon,periodic=False)

        lig_1 = t.topology.select('name C1 and resid 319')[0]
        lig_1_pair = [[lig_1,t.topology.select('name CA and resid ' + str(residue-1))[0]] for residue in binding_pocket_residue]

        lig_1_O = t.topology.select('name O1 and resid 319')[0]
        lig_1_O_pair = [[lig_1_O,t.topology.select('name CA and resid ' + str(residue-1))[0]] for residue in binding_pocket_residue]

        lig_1_com = t.topology.select('(type C or type O or type N) and resid 319')
        lig_1_com_cor = np.mean(t.xyz[:,lig_1_com,:],axis=1)

        lig_2 = t.topology.select('name C1 and resid 320')[0]
        lig_2_pair = [[lig_2,t.topology.select('name CA and resid ' + str(residue-1))[0]] for residue in binding_pocket_residue]

        lig_2_O = t.topology.select('name O1 and resid 320')
        lig_2_O_pair = [[lig_2_O,t.topology.select('name CA and resid ' + str(residue-1))[0]] for residue in binding_pocket_residue]

        lig_2_com = t.topology.select('(type C or type O or type N) and resid 320')
        lig_2_com_cor = np.mean(t.xyz[:,lig_2_com,:],axis=1)

        dist[:,(len(pairs)):(len(pairs)+len(lig_1_pair))] = md.compute_distances(t,lig_1_pair,periodic=False)

        dist[:,(len(pairs)+len(lig_1_pair)):(len(pairs)+2*len(lig_1_pair))] = md.compute_distances(t,lig_1_O_pair,periodic=False)

        com_pair_1 = []
        for j in range((len(pairs)+2*len(lig_1_pair)),(len(pairs)+3*len(lig_1_pair))):
            com_pair_1.append(['lig_1',t.topology.select('name CA and resid ' + str(binding_pocket_residue[j-len(pairs)-2*len(lig_1_pair)]))[0]])
            cen_mas  = t.xyz[:,t.topology.select('name CA and resid ' + str(binding_pocket_residue[j-len(pairs)-2*len(lig_1_pair)]))[0],:]
            dist[:,j] = np.linalg.norm((cen_mas-lig_1_com_cor),ord=2,axis=1)

        dist[:,(len(pairs)+3*len(lig_1_pair)):(len(pairs)+4*len(lig_1_pair))] = md.compute_distances(t,lig_2_pair,periodic=False)

        dist[:,(len(pairs)+4*len(lig_1_pair)):(len(pairs)+5*len(lig_1_pair))] = md.compute_distances(t,lig_2_O_pair,periodic=False)

        com_pair_2 = []

        for j in range((len(pairs)+5*len(lig_1_pair)),(len(pairs)+6*len(lig_1_pair))):
            com_pair_2.append(['lig_2',t.topology.select('name CA and resid ' + str(binding_pocket_residue[j-len(pairs)-5*len(lig_1_pair)]))[0]])
            cen_mas  = t.xyz[:,t.topology.select('name CA and resid ' + str(binding_pocket_residue[j-len(pairs)-5*len(lig_1_pair)]))[0],:]
            dist[:,j] = np.linalg.norm((cen_mas-lig_2_com_cor),ord=2,axis=1)

        newfile = filename.split('-strip.dcd')[0].split('./striped/')[1]
        p
        np.save('./numpy/'+newfile+'_msm_feature',dist)

        l = alpha_carbon + lig_1_pair + lig_1_O_pair + com_pair_1 + lig_2_pair + lig_2_O_pair + com_pair_2

        pickle.dump(l,open('msm_atom_pair.pkl','wb'))
        sys.exit()
    """
        #except:
        #    print(filename)
        #    os.system('cp ' + filename + ' ./striped_not_working/')
