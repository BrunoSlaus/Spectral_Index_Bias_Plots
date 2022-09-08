# -*- coding: utf-8 -*-
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Created on Sun Jul 16 2017

DESCRIPTION: This code plots the alpha (spectral index) vs the flux
(Total_flux_corr). This plot can be used to illustrate the bias
arising from the different detection limits of the two surveys used
for the calculation of the spectral indices. Since the NVSS survey 
has a higher limit, we add it to the plot.

The visual aspects of the plot should be altered dirrectly in the 
code: in the "plt.figure" part of the code...

Code Author: Bruno Slaus
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
##############################################################################################
#   User defined parameters:                                                                 #
##############################################################################################
basename  = '/home/bruno/Desktop/DIPLOMSKI_2/1_BIAS_Alpha/'    #Directory of the used catalogs
Det_Limit = 0.0025                                          #Detection limit of the NVSS field 

##############################################################################################
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits



#First we open the required data with astropy
print 'Plotting both the inner and the outer part of the field.\n'                                

Data=fits.open(basename + 'OUT_Single.fits')[1].data   
Total_flux_corr = Data['Total_flux_corr']
Spectral_index  = Data['Spectral_index']
print  'len(Total_flux_corr):', len(Total_flux_corr)
print  'len(Spectral_index):' , len(Spectral_index)

Data_2=fits.open(basename + 'IN_Single.fits')[1].data   
Total_flux_corr_2 = Data_2['Total_flux_corr']
Spectral_index_2  = Data_2['Spectral_index']
print  'len(Total_flux_corr_2):', len(Total_flux_corr_2)
print  'len(Spectral_index_2):' , len(Spectral_index_2)



#We define the limit line
x = np.arange(0.00001, 1, 0.0001)
alpha_limit = -(np.log10(x)-np.log10(Det_Limit)) / (np.log10(610)-np.log10(1400))



#Plotting the graph
fig = plt.figure(figsize=(7,7))
ax = fig.add_axes([0.12, 0.10, 0.8, 0.40])
ax.set_xlabel(r'$S_{610}/\ \mathrm{mJy}$', size=22)
ax.set_ylabel(r'$\alpha$', size=22)
ax.set_xlim(2*10**(-4)*1000, 10*1000)
ax.set_ylim(-2.5, 2.5)
plt.scatter(Total_flux_corr*1000, Spectral_index,color='black', s=2)
plt.plot(x*1000,alpha_limit, color='black', linestyle='--')
plt.xscale('log')
plt.hlines(0.745486, 0.0001, 10*1000, linestyle='-', color='red')
print 'Mean Outer Part: 0.745486'
#plt.hlines(0.87651, 0.0001, 10, linestyle=':', color='red')
#PLOTTING THE VERTICAL LINE OF THE LIMIT WE CUT/CROP
plt.vlines(0.02*1000, -3, 4, linestyle='-', color='blue')
ax.text(0.6*1000, 2, 'Vanjski dio', fontsize=16)

axx = fig.add_axes([0.12, 0.50, 0.8, 0.40])
#axx.set_xlabel(r'$S_{610}/\ mJy$', size=22)
axx.set_ylabel(r'$\alpha$', size=22)
axx.set_xlim(2*10**(-4)*1000, 10*1000)
axx.set_ylim(-1.9, 2)
plt.scatter(Total_flux_corr_2*1000, Spectral_index_2,color='black', s=2)
plt.plot(x*1000,alpha_limit, color='black', linestyle='--')
plt.xscale('log')
plt.hlines(0.646233, 0.0001, 10*1000, linestyle='-', color='red')
print 'Mean Inner Part: 0.646233'
#plt.hlines(0.81143, 0.0001, 10, linestyle=':', color='red')
#PLOTTING THE VERTICAL LINE OF THE LIMIT WE CUT/CROP
plt.vlines(0.015*1000, -3, 4, linestyle='-', color='blue')
axx.text(0.6*1000, 1.6, 'Srednji dio', fontsize=16)
axx.set_xticklabels([], alpha=0.0)

plt.savefig('Plot_Final', dpi=130)
#plt.savefig('/home/bruno/Desktop/BIAS_SEMINAR_TOGETHER', dpi=200)
plt.show()
plt.close()


#KRAJ

"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
End of code.
Modification history:


━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
