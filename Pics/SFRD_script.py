import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ax
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator
from matplotlib import rc
rc('font',**{'family':'serif'})

#specify file stuff here
filename_raw = "SFRD_raw.txt"
filename_src = "SFRD_src.txt"
file_raw = open(filename_raw, 'r')
file_src = open(filename_src, 'r')

#declare main arrays
z_raw = np.array([])
sfrd_raw = np.array([])
z_src = np.array([])
sfrd_src = np.array([])

for line in file_raw.readlines():
    cols = line.split()
    z_raw = np.append(z_raw, float(cols[0]))
    sfrd_raw = np.append(sfrd_raw, float(cols[1]))

for line in file_src.readlines():
    cols = line.split()
    z_src = np.append(z_src, float(cols[0]))
    sfrd_src = np.append(sfrd_src, float(cols[1]))

#plot the two vs each other
fig,ax1 = plt.subplots(figsize=(12, 12))
plt.plot(z_raw, sfrd_raw, color='r', linewidth=5, label='Park et al.')
plt.plot(z_src, sfrd_src, color='b', linewidth=5, label='LIMFAST')
plt.xlabel('$\\rm redshift$', fontsize=30)
plt.ylabel('${\\rm SFRD\\,\\,{ [M_{\\odot}\\,yr^{-1}\\,Mpc^{-3}]}}$', fontsize=30)
ax = plt.gca()
ax1.set_xlim(6.0, 16.0)
ax1.set_ylim(10e-6, 1)
ax1.set_xscale('linear')
ax1.set_yscale('log')
ax1.xaxis.set_major_locator(MultipleLocator(2))
ax1.xaxis.set_minor_locator(MultipleLocator(0.5))
ax1.tick_params(axis='y',which='major',labelsize=25,length=14,width=5.)
ax1.tick_params(axis='y',which='minor',labelsize=25,length=11,width=2)
ax1.tick_params(axis='x',which='major',labelsize=25,length=14,width=5.)
ax1.tick_params(axis='x',which='minor',labelsize=25,length=11,width=2)
plt.legend(fontsize=28,frameon=False,loc='best')
plt.savefig("SFRD.png",bbox_inches='tight')
