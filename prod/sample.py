import glob
import numpy as np
import mdtraj as md
import sys
top = '../CB1_assym-stripped.prmtop'
totdist = []
totfile = []
totframes = []
for filename in glob.glob('./striped/CB1_assym_pr_8_*'):
  t = md.load(filename,top=top)
  nframes = t.n_frames
  #ASP = t.topology.select('name OG and resid 294')
  resid = t.topology.select('name CA and ((resid 288 to 298) or (resid 27 to 37))')
  cen_mas = np.mean(t.xyz[:,resid,:],axis=1)
  dist = np.empty([nframes,2])
  LIG = t.topology.select('name C1 and resname LIG')
  #dist = md.compute_distances(t,[[ASP[0],LIG[0]],[ASP[0],LIG[1]]],periodic=False)
  dist[:,0] = np.linalg.norm(t.xyz[:,LIG[0],:]-cen_mas,ord=2,axis=1)
  #print(dist.shape,t.xyz[:,LIG[1],:].shape,cen_mas.shape)
  dist[:,1] = np.linalg.norm(t.xyz[:,LIG[1],:]-cen_mas,ord=2,axis=1)
  #sys.exit()
  #print(np.min(dist,axis=1))
  files = [filename for i in range(len(dist))]
  frames = [i for i in range(len(dist))]
  dist = np.min(dist,axis=1)
  totdist.append(dist)
  totfile.append(files)
  totframes.append(frames)
txx = list(np.concatenate(totdist))
txxfile = list(np.concatenate(totfile))
txxframes = list(np.concatenate(totframes))
a,files,frames = zip(*sorted(zip(txx,txxfile,txxframes)))
print(a[:15])
for i,file in enumerate(files[:10]):
  f = open('../cpp/CB1_assym_pr_9_frame_'+str(i),'w')
  newfile = file.replace('-strip.dcd','.nc').split('./striped/')[1]
  f.write('parm ../CB1_assym.prmtop \n')
  f.write('trajin ../traj/' +newfile + '\n')
  f.write('autoimage origin \n')
  f.write('center origin \n')
  f.write('trajout ' + '../rst/CB1_assym_pr_8_frame_'+str(i+26) + '.rst7 onlyframes ' + str(frames[i]+1) + '\n')
  f.write('run \n')
  f.write('quit \n')
  f.close()
