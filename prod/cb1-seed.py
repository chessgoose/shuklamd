
import numpy as np
import mdtraj as md
import glob
import sys
import random
tot_trp_z = []
tot_lig_z = []
tot_tyr_z = []
totfile = []
totdist = []
zcoords = []

binding_pocket_residue = [170,174,177,193,196,197,267,268,269,271,275,279,359,363,379,383]
for filename in glob.glob('./numpy/CB1_assym_pr_8_frame_*_lig_1.npy'):
    #t = md.load(filename,top='CB2_assym-strip.prmtop')
    dist = np.load(filename)
    totdist.append(dist[:,0])
    zcoords.append(dist[:,4])
    totfile.append(filename)

traj = {}
required_frame = 10
txx =  np.concatenate(totdist)

print(sorted(txx)[:15])
x_min = np.min(txx)
#sys.exit()
for file,x_data,y_data in zip(totfile,totdist,zcoords):
    x_frame = set(np.where(x_data>1)[0].tolist()) & set(np.where(x_data<1.5)[0].tolist())
    y_frame = set(np.where(y_data>0)[0].tolist()) & set(np.where(y_data<0.5)[0].tolist())
    #z_frame = set(np.where(z_data>(z_max/2))[0].tolist()) & set(np.where(z_data<z_max)[0].tolist())
    frame = list(x_frame & y_frame)
    print(frame,file)
    f = file.split('CB1_assym_pr_')[1]
    fi = f.split('_frame_')[0]
    f_new = f.split('_frame_')[1].split('_lig_1.npy')[0]
    if len(frame) > 0:
        #print(fi,f_new,frame)
        keys = list(traj.keys())
        if fi in keys:
            traj[fi].update({f_new:frame})
        else:
            traj[fi] = {f_new:frame}
#sys.exit()
for i in range(0,required_frame):
     keys = list(traj.keys())
     #print(keys)
     if keys == []:
         break
     random.seed()
     random.shuffle(keys)
     indi = list(traj[keys[0]].keys())
     random.shuffle(indi)
     frames = traj[keys[0]][indi[0]]
     #print(frames)
     random.shuffle(frames)
     frame = frames[0]
     traj[keys[0]][indi[0]].remove(frame)
     if traj[keys[0]][indi[0]] == []:
         del traj[keys[0]][indi[0]]
     if not bool(traj[keys[0]]):
         del traj[keys[0]]
     #print(traj)
     f = open('../cpp/CB1_assym_pr_9_frame_' +str(i)+'_cpp','w')
     f.write('parm ../CB1_assym.prmtop \n')
     f.write('trajin ../traj/CB1_assym_pr_'+ keys[0] + '_frame_' + indi[0] + '.nc \n')
     f.write('trajout ../rst/CB1_assym_pr_9_' +str(i)+'.rst7 onlyframes '+ str(frame +1) + ' \n')
     f.write('run \n')
     f.write('quit \n')
