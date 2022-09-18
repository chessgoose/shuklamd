
import numpy as np
import mdtraj as md
import glob
import sys
import random
import pyemma

totfile = []
totdata = []
f_num = 500

for filename in glob.glob('./numpy/c_CB2_assym_pr_*_frame_*_lig_1.npy'):
    dist = np.load(filename)
    f = filename.split('c_')[1]
    totdata.append(dist) #append ndarray
    totfile.append(f.split('_lig_1.npy')[0])

#Append all trajectory arrays, input into clustering
cluster_k = pyemma.coordinates.cluster_kmeans(totdata,k=f_num,max_iter=1000,stride=5)
o_frames = cluster_k.get_output()

required_frames = []
roun = "1"
number = 1

for i in range(0, f_num):
    z = []
    for j,dj in enumerate(o_frames):
        if i in dj:
           a = np.where(dj==i)
           temp  = [[j,ind] for ind in a[0]]
           z.extend(temp)
    random_permu = np.random.permutation(len(z))
    frame = [[totfile[z[j][0]],z[j][1]] for j in random_permu[:number]]
    required_frames.extend(frame)

for i in range(len(required_frames)):
    f = open('../cpp/c_CB2_assym_pr_'+roun+'_frame_'+str(i),'w')
    f.write('parm ../CB2_assym.prmtop \n')
    f.write('trajin ../traj/' + required_frames[i][0] + '.nc \n')
    f.write('trajout ../rst/'+ 'c_CB2_assym_pr_' + roun+ '_frame_' + str(i+1) + '.rst7 onlyframes ' + str(required_frames[i][1]+1) + ' \n')
    f.write('run \n')
    f.write('quit \n')
    f.close()
