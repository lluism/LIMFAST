import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ax
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator
from matplotlib import rc
rc('font',**{'family':'serif'})

#file = open('./path_to_file_power_spectrum/name_file','r')
filename1 = "../Pics/LineLum_HI6563A_z012.00_PS_Mom"
filename2 = "../Pics/LineLum_HI6563A_z012.00_PS_Ene"
filename3 = "../Pics/LineLum_HI6563A_z005.00_PS_Mom"
filename4 = "../Pics/LineLum_HI6563A_z005.00_PS_Ene"

#filename1 = "../Pics/LineLum_OIII5007A_z010.00_PS_Mom"
#filename2 = "../Pics/LineLum_OIII5007A_z010.00_PS_Ene"
#filename3 = "../Pics/LineLum_OIII5007A_z006.00_PS_Mom"
#filename4 = "../Pics/LineLum_OIII5007A_z006.00_PS_Ene"

k1 = np.array([])
k2 = np.array([])
k3 = np.array([])
k4 = np.array([])

pk1 = np.array([])
pk2 = np.array([])
pk3 = np.array([])
pk4 = np.array([])

delta_k1 = np.array([])
delta_k2 = np.array([])
delta_k3 = np.array([])
delta_k4 = np.array([])

file1 = open(filename1)
for lines1 in file1.readlines():
    cols1 = lines1.split()
    k1 = np.append(k1, float(cols1[0]))
    pk1 = np.append(pk1, float(cols1[1]))
    delta_k1 = np.append(delta_k1, float(cols1[2]))

file2 = open(filename2)
for lines2 in file2.readlines():
    cols2 = lines2.split()
    k2 = np.append(k2, float(cols2[0]))
    pk2 = np.append(pk2, float(cols2[1]))
    delta_k2 = np.append(delta_k2, float(cols2[2]))

file3 = open(filename3)
for lines3 in file3.readlines():
    cols3 = lines3.split()
    k3 = np.append(k3, float(cols3[0]))
    pk3 = np.append(pk3, float(cols3[1]))
    delta_k3 = np.append(delta_k3, float(cols3[2]))

file4 = open(filename4)
for lines4 in file4.readlines():
    cols4 = lines4.split()
    k4 = np.append(k4, float(cols4[0]))
    pk4 = np.append(pk4, float(cols4[1]))
    delta_k4 = np.append(delta_k4, float(cols4[2]))

# PLOT DELTA(k)
fig,ax1 = plt.subplots(figsize=(12,12))
plt.plot(k1,delta_k1/(1.0e-23)**2,ls='-.',color='k',linewidth=2, label='Momentum')
plt.plot(k2,delta_k2/(1.0e-23)**2,ls='-',color='k',linewidth=2, label='Energy')
plt.plot(k3,delta_k3/(1.0e-23)**2,ls='-.',color='g',linewidth=2)
plt.plot(k4,delta_k4/(1.0e-23)**2,ls='-',color='g',linewidth=2)
HALPHAPS_GONG_Z5 = np.loadtxt('HaPS_Gong_z5.csv', delimiter=',')
plt.plot(HALPHAPS_GONG_Z5[:,0], HALPHAPS_GONG_Z5[:,1], 'r:', label=r'$\Delta_\mathrm{H\alpha}^2(k,z=5)$ from Gong+17')
plt.plot(HALPHAPS_GONG_Z5[:,0], HALPHAPS_GONG_Z5[:,1]*100., 'r:', alpha=0.3)
plt.xlabel('$k\ {\\rm [h\\,  Mpc^{-1}]}$',fontsize=30);
# units below need to be corrected
plt.ylabel('${\\Delta_\mathrm{H\\alpha}^2 (k)\ {\\rm [(Jy/sr)^2]}}$',fontsize=30)
ax = plt.gca()
#ax1.set_xlim(0.,22.)
ax1.set_yscale('log')
ax1.set_xscale('log')
#ax1.xaxis.set_major_locator(MultipleLocator(10))
#ax1.xaxis.set_minor_locator(MultipleLocator(5))
ax1.yaxis.set_major_locator(LogLocator(base = 10.0))
ax1.xaxis.set_major_locator(LogLocator(base = 10.0))
majorFormatter = FormatStrFormatter('%d')
ax1.tick_params(axis='y',which='major',labelsize=25,length=14,width=5.)
ax1.tick_params(axis='y',which='minor',labelsize=25,length=11,width=2)
ax1.tick_params(axis='x',which='major',labelsize=25,length=14,width=5.)
ax1.tick_params(axis='x',which='minor',labelsize=25,length=11,width=2)
plt.legend(fontsize=25,frameon=False,numpoints=1,loc='lower right')
#plt.show()
#plt.savefig(filename+'_Deltak.pdf', bbox_inches='tight')
plt.savefig('HI6563A_Feedback_Deltak.pdf', bbox_inches='tight')
#plt.savefig('OIII5007A_Feedback_Deltak.pdf', bbox_inches='tight')
