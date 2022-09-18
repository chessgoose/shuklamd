import numpy as np
import glob
import pyemma
import pickle
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib as mpl

totdist = []
filename = []
roun = [i for i in range(0,2000)]
for i in roun:
    for file in glob.glob('./numpy/CB2_assym_pr_'+str(i)+'_frame_*_lig_com.npy'):
        distI = np.load(file)
        totdist.append(distI)
        #filI = file.split('./numpy/')[1]
        #filI = filI.split('_rrcs.npy')[0]
        #filename.append('./CB2-APO_inactive/pmemd_testing/analysis/striped/'+filI + '.nc')
        del distI

txx = np.concatenate(totdist)
print(txx.shape)
#txx57 = np.concatenate(totTYR57)
#parametization
x_bins = 150
y_bins = 150
R = 0.001987
T = 300
fig_wid = 10
fig_hig = 7
cmap = mpl.cm.jet
Max_energy =  6
flag=0
weights=None

x_data =  txx[:,j]*10
print(min(x_data))
#x_data = txx57[:,0]*10
y_data =  txx[:,j+1]*10
x_data_min =  np.min(x_data)
y_data_min =  np.min(y_data)
x_data_max =  np.max(x_data)
y_data_max =  np.max(y_data)
x_hist_lim_low =  x_data_min -0.5
y_hist_lim_low =  y_data_min -0.5
x_hist_lim_high = x_data_max +0.5
y_hist_lim_high = y_data_max  +0.5
x_lim_low = (int(np.min(x_data)/5))*5.0
y_lim_low = (int(np.min(y_data)/5))*5.0
x_lim_high = (int(np.max(x_data)/5))*5.0
y_lim_high = (int(np.max(y_data)/5)+1)*5.0
hist= np.histogram2d(x_data,y_data, bins=[x_bins,y_bins],
                     range = [[x_hist_lim_low,x_hist_lim_high],[y_hist_lim_low,y_hist_lim_high]],
                     density= True,weights=weights)
prob_density = hist[0]
xedge = hist[1]
yedge = hist[2]
x_bin_size = xedge[1]-xedge[0]
y_bin_size = yedge[1]-yedge[0]
free_energy = -R*T*np.log(prob_density*x_bin_size*y_bin_size)
min_free_energy= np.min(free_energy)
delta_free_energy = free_energy - min_free_energy
xx = [(xedge[i]+xedge[i+1])/2 for i in range(len(xedge)-1)]
yy = [(yedge[i]+yedge[i+1])/2 for i in range(len(yedge)-1)]
fig, axs = plt.subplots(1,1,figsize=(fig_wid,fig_hig))
cd =axs.contourf(xx,yy,delta_free_energy.T, np.linspace(0,Max_energy,Max_energy*5+1),
                 vmin=0.0, vmax=Max_energy,cmap=cmap)
cbar = fig.colorbar(cd,ticks=range(Max_energy+1))
cbar.ax.set_yticklabels(range(Max_energy+1),fontsize=22)
cbar.ax.set_ylabel("Free energy", labelpad=15,**hfont,fontsize=25)
# axs.scatter(inactive[0,dica[x_key]]*10,inactive[0,dica[y_key]]*10,s=100,c='blue',marker='^')
# axs.scatter(active[0,dica[x_key]]*10,active[0,dica[y_key]]*10,s=100,c='red',marker='o')
axs.set_xlim([x_lim_low,x_lim_high])
axs.set_ylim([y_lim_low,y_lim_high])
axs.set_xlim([-40,40])
axs.set_ylim([0,80])
#axs.set_xticks(range(int(x_lim_low),int(x_lim_high)+1,5))
#axs.set_xticklabels(range(int(x_lim_low),int(x_lim_high)+1,5))
axs.set_xticks(range(int(-40),int(40)+1,20))
axs.set_xticklabels(range(int(-40),int(40)+1,20))
axs.set_yticks(range(int(0),int(80)+1,20))
axs.set_yticklabels(range(int(0),int(80)+1,20))
#axs.set_yticks(range(int(5),int(20)+1,5))
#axs.set_yticklabels(range(int(5),int(20)+1,5))
plt.xlabel('x lig1 (\AA)', **hfont,fontsize=30)
plt.ylabel('z lig1 (\AA)', **hfont,fontsize=30)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.tight_layout()
#plt.grid()
plt.savefig('lig_binding_lig1_com_x',transparent=True,dpi=500)
