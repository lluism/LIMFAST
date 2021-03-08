import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ax
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator
from matplotlib import rc
rc('font',**{'family':'serif'})

#file = open('./path_to_file_power_spectrum/name_file','r')
filename = "../Out-of-box_Output/ps_Lyaigm_z007.08_HIIfilter1_RHIImax50_200_300Mpc"
file = open(filename)
k = np.array([])
pk = np.array([])
delta_k = np.array([])
for lines in file.readlines():
    cols = lines.split()
    k = np.append(k, float(cols[0]))
    pk = np.append(pk, float(cols[1]))
    delta_k = np.append(delta_k, float(cols[2]))

# PLOT P(K)

fig,ax1 = plt.subplots(figsize=(12,12))
plt.plot(k,pk,linestyle='solid',color='k',linewidth=6, label='P(k)')
plt.xlabel('$k {\\rm [h\\,  Mpc^{-1}]}$',fontsize=30);
# units below need to be corrected
plt.ylabel('${P(k) \\, {\\rm [erg\\,s^{-1}\\,cm^{-2}\\,arcsec^{-2}]}}$',fontsize=30)
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
plt.legend(fontsize = 28,frameon=False,numpoints=1,loc='best')
#plt.show()
plt.savefig(filename+'_k.pdf', bbox_inches='tight')

# PLOT DELTA(k)
fig,ax1 = plt.subplots(figsize=(12,12))
plt.plot(k,delta_k,linestyle='solid',color='k',linewidth=6, label='Delta_k')
plt.xlabel('$k {\\rm [h\\,  Mpc^{-1}]}$',fontsize=30);
# units below need to be corrected
plt.ylabel('${\\Delta^2 (k) \\, {\\rm [erg\\,s^{-1}\\,cm^{-2}\\,arcsec^{-2}]}}$',fontsize=30)
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
plt.legend(fontsize = 28,frameon=False,numpoints=1,loc='best')
#plt.show()
plt.savefig(filename+'_Deltak.pdf', bbox_inches='tight')
