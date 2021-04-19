import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ax
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator
from matplotlib import rc
rc('font',**{'family':'serif'})

#file = open('./path_to_file_power_spectrum/name_file','r')
#filename = "../Out-of-box_Output/ps_Lyaigm_z007.08_HIIfilter1_RHIImax50_200_300Mpc"
#filename1 = "../Pics/LineLum_HI6563A_z015.00_PS"
filename1 = "../Pics/LineLum_HI6563A_z005.00_PS"
k1 = np.array([])
pk1 = np.array([])
delta_k1 = np.array([])



#filename2 = "../Pics/LineLum_HI6563A_z010.00_PS_Ene"
'''
filename2 = "../Pics/LineLum_OIII5007A_z015.00_PS_Mom"
filename3 = "../Pics/LineLum_HI6563A_z005.00_PS_Mom"
filename4 = "../Pics/LineLum_OIII5007A_z005.00_PS_Mom"
filename5 = "../Pics/LineLum_HI6563A_z010.00_PS_Mom"
filename6 = "../Pics/LineLum_OIII5007A_z010.00_PS_Mom"

k1 = np.array([])
k2 = np.array([])
k3 = np.array([])
k4 = np.array([])
k5 = np.array([])
k6 = np.array([])

pk1 = np.array([])
pk2 = np.array([])
pk3 = np.array([])
pk4 = np.array([])
pk5 = np.array([])
pk6 = np.array([])

delta_k1 = np.array([])
delta_k2 = np.array([])
delta_k3 = np.array([])
delta_k4 = np.array([])
delta_k5 = np.array([])
delta_k6 = np.array([])
'''

file1 = open(filename1)
for lines1 in file1.readlines():
    cols1 = lines1.split()
    k1 = np.append(k1, float(cols1[0]))
    pk1 = np.append(pk1, float(cols1[1]))
    delta_k1 = np.append(delta_k1, float(cols1[2]))

'''
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

file5 = open(filename5)
for lines5 in file5.readlines():
    cols5 = lines5.split()
    k5 = np.append(k5, float(cols5[0]))
    pk5 = np.append(pk5, float(cols5[1]))
    delta_k5 = np.append(delta_k5, float(cols5[2]))

file6 = open(filename6)
for lines6 in file6.readlines():
    cols6 = lines6.split()
    k6 = np.append(k6, float(cols6[0]))
    pk6 = np.append(pk6, float(cols6[1]))
    delta_k6 = np.append(delta_k6, float(cols6[2]))
'''

# PLOT P(K)


fig,ax1 = plt.subplots(figsize=(12,12))
plt.plot(k1,pk1/(1.0e-23)**2,'ko-',linewidth=2, label='Ha: z=6')
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
plt.savefig('HI6563A_DEFAULTBOX_k.pdf', bbox_inches='tight')
plt.clf()


HALPHAPS_GONG_Z5 = np.loadtxt('HaPS_Gong_z5.csv', delimiter=',')

# PLOT DELTA(k)
fig,ax1 = plt.subplots(figsize=(12,12))
plt.plot(k1,delta_k1/(1.0e-23)**2,ls='--',color='k',linewidth=2, label='Ha: z=5')
#Pshot_approx1 = delta_k1[-1]/(1.0e-23)**2 / k1[-1]**3 * 2.*np.pi**2
#plt.plot(k1,Pshot_approx1*k1**3/2/np.pi**2,ls=':',color='k',linewidth=2)
#Pshot_approx2 = delta_k2[-1]/(1.0e-23)**2 / k2[-1]**3 * 2.*np.pi**2
#plt.plot(k2,Pshot_approx2*k2**3/2/np.pi**2,ls=':',color='b',linewidth=2)
'''
plt.plot(k1,delta_k1/(1.0e-23)**2/50.,ls=':',color='k',linewidth=1, label='Ha (rescaled)')
plt.plot(k2,delta_k2/(1.0e-23)**2,ls='-',color='k',linewidth=2, label='OIII')
plt.plot(k3,delta_k3/(1.0e-23)**2,ls='--',color='g',linewidth=2)
plt.plot(k3,delta_k3/(1.0e-23)**2/10.,ls=':',color='g',linewidth=1)
plt.plot(k4,delta_k4/(1.0e-23)**2,ls='-',color='g',linewidth=2)
plt.plot(k5,delta_k5/(1.0e-23)**2,ls='--',color='r',linewidth=2)
plt.plot(k5,delta_k5/(1.0e-23)**2/25.,ls=':',color='r',linewidth=1)
plt.plot(k6,delta_k6/(1.0e-23)**2,ls='-',color='r',linewidth=2)
plt.text(8.0e-1, 4.0e-4, r'$z=15$', color='k', fontsize=18)
plt.text(8.0e-1, 2.0e-4, r'$z=10$', color='r', fontsize=18)
plt.text(8.0e-1, 1.0e-4, r'$z=5$', color='g', fontsize=18)
'''
plt.plot(HALPHAPS_GONG_Z5[:,0], HALPHAPS_GONG_Z5[:,1], 'r:', label=r'$\Delta_\mathrm{H\alpha}^2(k,z=5)$ from Gong+17')
plt.plot(HALPHAPS_GONG_Z5[:,0], HALPHAPS_GONG_Z5[:,1]*100., 'r:', label=r'$\Delta_\mathrm{H\alpha}^2(k,z=5)$ from Gong+17')
plt.xlabel('$k\ {\\rm [h\\,  Mpc^{-1}]}$',fontsize=30);
# units below need to be corrected
plt.ylabel('${\\Delta_\mathrm{H\\alpha}^2 (k)\ {\\rm [(Jy/sr)^2]}}$',fontsize=30)
ax = plt.gca()
#ax1.set_xlim(0.,22.)
ax1.set_yscale('log')
ax1.set_xscale('log')
#ax1.set_ylim([3.0e-4, 6.0e1])
#ax1.xaxis.set_major_locator(MultipleLocator(10))
#ax1.xaxis.set_minor_locator(MultipleLocator(5))
ax1.yaxis.set_major_locator(LogLocator(base = 10.0))
ax1.xaxis.set_major_locator(LogLocator(base = 10.0))
majorFormatter = FormatStrFormatter('%d')
ax1.tick_params(axis='y',which='major',labelsize=25,length=14,width=5.)
ax1.tick_params(axis='y',which='minor',labelsize=25,length=11,width=2)
ax1.tick_params(axis='x',which='major',labelsize=25,length=14,width=5.)
ax1.tick_params(axis='x',which='minor',labelsize=25,length=11,width=2)
plt.legend(fontsize=25,frameon=False,numpoints=1,loc='upper left')
#plt.show()
#plt.savefig(filename+'_Deltak.pdf', bbox_inches='tight')
#plt.savefig('HI6563A_OIII5007A_Deltak.pdf', bbox_inches='tight')
plt.savefig('HI6563A_DEFAULTBOX_Deltak.pdf', bbox_inches='tight')
