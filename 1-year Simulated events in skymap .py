#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import astropy.units as u
from astropy.cosmology import Planck18 as cosmo
import pandas as pd
from scipy import stats
from firesong.Firesong import firesong_simulation
from firesong.Evolution import get_evolution, SourcePopulation
from firesong.Evolution import TransientSourcePopulation, cosmology
from firesong.Luminosity import get_LuminosityFunction
import csv
from astropy.coordinates import SkyCoord
from astropy.visualization import astropy_mpl_style
from tabulate import tabulate


# In[2]:


#luminosity:erg/yr, Luminosity(Optic), luminosizy(neutrinos)
#Timescale = 100-1000
#E_tot = 2e45 #erg
#eta = 2.7e3
Timescale = 1000 #s
#Luminosity = E_tot/Timescale #optic,erg/S
#L_nu = Luminosity*365*24*60*60*eta #neutrinos,erg/yr
L_nu = 1.7e53
zmax = 2
index=2.28
density = 2.58e-5 #local density https://doi.org/10.1111/j.1365-2966.2011.18162.x 
#local volumetric rate of 0.301 ± 0.062, 0.258 ± 0.072 and 0.447 ± 0.139 for SNe Ia, Ibc and II,respectively (in units of 10−4 SN Mpc−3 yr−1).  


sim = firesong_simulation(None,filename=None, 
                          density = density,                           
                          Transient = 'Ture', #just take Ibc
                          Evolution = 'MD2014SFR',
                          timescale = Timescale, #100-1000
                          luminosity = L_nu, #erg/yr
                          emin=1e4,
                          emax=1e7,
                          index=index,
                          zmax = zmax)


# In[3]:


"""
Data from FIRESONG
sim['sources']
dec : Random declination over the entire sky (deg)
ra  : Right ascension (deg)
flux: (GeV/cm^2.s.sr)
z   : redshift
lumis : (erg/yr)
D_L: Luminosity distance(Mpc)
"""


flux = sim['sources']['flux']#neutrino flux
dec = sim['sources']['dec']

days = 365 # /yr

#Total diffuse flux
#fluence is the flux without transient mode
#flux is the FIRESONG output

zs = sim['sources']['z']
population = TransientSourcePopulation(cosmology,get_evolution("MD2014SFR"),timescale=Timescale)
D_L = population.LuminosityDistance(zs)  #Luminosity distance in Mpc
fluence = flux*((1.+zs)*1000)
TOT = sum(fluence)
TOT = TOT/population.yr2sec #/s
TOT /= 4*np.pi # /sr
print('total diffuse flux',TOT,'Gev/cm^2.s.sr')


# In[4]:


# Append A_eff array to the list
gfile = np.genfromtxt('Effa_all_streams_gold_only.txt')
genergy_bins = gfile[:,0]
g_eff_list = []  # List to store A_eff arrays

for j in range(len(dec)):
    if dec[j] > 30:
        A_eff = 1
    elif 30 >= dec[j] > 0:
        A_eff = 2
    elif 0 >= dec[j] > -5:
        A_eff = 3
    elif -5 >= dec[j] > -30:
        A_eff = 4
    else:
        A_eff = 5

    g_eff_list.append(A_eff)

j=np.array(g_eff_list)

gnus = 0
for i in range(0,len(genergy_bins)): 

    gnus += gfile[:,j][i] *1e4*(10**5)**(index-2)*fluence /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print('expected number of events of signel sourse:',gnus)
print('expected number of events:',sum(gnus))
print('maximal of signel sourse:',max(gnus),
'minimal:',min(gnus))


# In[5]:


#cautulate the nus from totol flux

#dec 30-90 
g_eff_fttfl_1 = gfile[:, 1] 
gFlux_fttfl_1 = 2*np.pi*TOT*(np.cos(np.pi*2/3)-np.cos(np.pi))

#dec 0-30
g_eff_fttfl_2 = gfile[:, 2] 
gFlux_fttfl_2 = 2*np.pi*TOT*(np.cos(np.pi*1/2)-np.cos(np.pi*2/3))

#dec -5-0
g_eff_fttfl_3 = gfile[:, 3]
gFlux_fttfl_3 = 2*np.pi*TOT*(np.cos(np.pi*17/36)-np.cos(np.pi*1/2))

#dec -30--5
g_eff_fttfl_4 = gfile[:, 4] 
gFlux_fttfl_4 = 2*np.pi*TOT*(np.cos(np.pi*1/3)-np.cos(np.pi*17/36))

#dec -90--30
g_eff_fttfl_5 = gfile[:, 5]
gFlux_fttfl_5 = 2*np.pi*TOT*(np.cos(np.pi*0)-np.cos(np.pi*1/3))

gnus_fttfl_1 = 0
for i in range(0,len(genergy_bins)):  
    gnus_fttfl_1 += g_eff_fttfl_1[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_1 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(gnus_fttfl_1,'dec 30-90')
gnus_fttfl_2 = 0
for i in range(0,len(genergy_bins)):  
    gnus_fttfl_2 += g_eff_fttfl_2[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_2 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(gnus_fttfl_2,'dec 30-90')
gnus_fttfl_3 = 0
for i in range(0,len(genergy_bins)):  
    gnus_fttfl_3 += g_eff_fttfl_3[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_3 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(gnus_fttfl_3,'dec 30-90')
gnus_fttfl_4 = 0
for i in range(0,len(genergy_bins)):  
    gnus_fttfl_4 += g_eff_fttfl_4[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_4 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(gnus_fttfl_4,'dec 30-90')
gnus_fttfl_5 = 0
for i in range(0,len(genergy_bins)):  
    gnus_fttfl_5 += g_eff_fttfl_5[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_5 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(gnus_fttfl_5,'dec 30-90')

print('expected number from total flux Gold: ' , gnus_fttfl_1+gnus_fttfl_2+gnus_fttfl_3+gnus_fttfl_4+gnus_fttfl_5,'/yr')


# In[6]:


# Append A_eff array to the list
afile = np.genfromtxt('Effa_all_streams_gold_bronze.txt')
genergy_bins = gfile[:,0]
bfile = afile-gfile

bnus=0
for i in range(0,len(genergy_bins)): 

    bnus += bfile[:,j][i] *1e4*(10**5)**(index-2)*fluence /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print('expected number of events of signel sourse bronze:',bnus)
print('expected number of events:',sum(bnus))
print('maximal of signel sourse:',max(bnus),
'minimal:',min(bnus))


# In[7]:


#cautulate the nus from totol flux

#dec 30-90 
b_eff_fttfl_1 = bfile[:, 1] 
bFlux_fttfl_1 = 2*np.pi*TOT*(np.cos(np.pi*2/3)-np.cos(np.pi))

#dec 0-30
b_eff_fttfl_2 = bfile[:, 2] 
bFlux_fttfl_2 = 2*np.pi*TOT*(np.cos(np.pi*1/2)-np.cos(np.pi*2/3))

#dec -5-0
b_eff_fttfl_3 = bfile[:, 3]
bFlux_fttfl_3 = 2*np.pi*TOT*(np.cos(np.pi*17/36)-np.cos(np.pi*1/2))

#dec -30--5
b_eff_fttfl_4 = bfile[:, 4] 
bFlux_fttfl_4 = 2*np.pi*TOT*(np.cos(np.pi*1/3)-np.cos(np.pi*17/36))

#dec -90--30
b_eff_fttfl_5 = bfile[:, 5]
bFlux_fttfl_5 = 2*np.pi*TOT*(np.cos(np.pi*0)-np.cos(np.pi*1/3))

bnus_fttfl_1 = 0
for i in range(0,len(genergy_bins)):  
    bnus_fttfl_1 += b_eff_fttfl_1[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_1 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(bnus_fttfl_1,'dec 30-90')
bnus_fttfl_2 = 0
for i in range(0,len(genergy_bins)):  
    bnus_fttfl_2 += b_eff_fttfl_2[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_2 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(bnus_fttfl_2,'dec 30-90')
bnus_fttfl_3 = 0
for i in range(0,len(genergy_bins)):  
    bnus_fttfl_3 += b_eff_fttfl_3[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_3 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(bnus_fttfl_3,'dec 30-90')
bnus_fttfl_4 = 0
for i in range(0,len(genergy_bins)):  
    bnus_fttfl_4 += b_eff_fttfl_4[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_4 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(bnus_fttfl_4,'dec 30-90')
bnus_fttfl_5 = 0
for i in range(0,len(genergy_bins)):  
    bnus_fttfl_5 += b_eff_fttfl_5[i]*1e4*86400 * days * (10**5)**(index-2)*gFlux_fttfl_5 /(1-index)*((genergy_bins[i])**(1-index)-(genergy_bins[i-1])**(1-index))
print(bnus_fttfl_5,'dec 30-90')

print('expected number from total flux Bronze: ' , bnus_fttfl_1+bnus_fttfl_2+bnus_fttfl_3+bnus_fttfl_4+bnus_fttfl_5,'/yr')


# In[8]:


# Assuming you have a probability array
g_prob_array = gnus

# Generate random Poisson-distributed values for each element in the probability array
g_lambda_values = -np.log(1 - g_prob_array)  # Use inverse transform sampling

# Generate outcomes based on Poisson-distributed values and a threshold
threshold = 1  # You can choose your own threshold
g_binary_array = (np.random.poisson(g_lambda_values) >= threshold).astype(int)
g_indices = np.where(g_binary_array == 1)[0]

print("Indices of Gold Alerts:", g_indices)


# In[9]:


threshold2 = 2  # You can choose your own threshold
g2_binary_array = (np.random.poisson(g_lambda_values) >= threshold2).astype(int)
g2_indices = np.where(g2_binary_array == 1)[0]

print("Indices of Gold Alerts:", g2_indices)


# In[10]:


# Assuming you have a probability array
b_prob_array = bnus

# Generate random Poisson-distributed values for each element in the probability array
b_lambda_values = -np.log(1 - b_prob_array)  # Use inverse transform sampling

# Generate outcomes based on Poisson-distributed values and a threshold
threshold = 1  # You can choose your own threshold
b_binary_array = (np.random.poisson(b_lambda_values) >= threshold).astype(int)
b_indices = np.where(b_binary_array == 1)[0]

print("Indices of Bronze Alerts:", b_indices)


# In[11]:


threshold2 = 2  # You can choose your own threshold
b2_binary_array = (np.random.poisson(b_lambda_values) >= threshold2).astype(int)
b2_indices = np.where(b2_binary_array == 1)[0]

print("Indices of bronze Alerts:", b2_indices)


# In[12]:


mean, std = -17.68288114746782, 1.0472514600441147
mu, sigma = mean, std # mean and standard deviation
s = np.random.normal(mu, sigma, int(population.Nsources(density, zmax)))
abs(mu - np.mean(s))
abs(sigma - np.std(s, ddof=1))
count, bins, ignored = plt.hist(s, 500, density=True)
plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
               np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
         linewidth=2, color='r')
plt.show()
Bol_MAG=s


# In[13]:


#Bol_MAG = -17.683136246786635 #assume absolute Mag standard candel
def mag_bol(Bol_MAG,D_L):
    mag_bol = Bol_MAG-5*np.log10(10/D_L*10**-6)
    return(mag_bol)
mag = mag_bol(Bol_MAG,D_L)


# In[14]:


bronze = (0.8, 0.5, 0.2) 
ra=sim['sources']['ra']-180
b_dec=dec[b_indices]
g_dec=dec[g_indices]
b_mag=mag[b_indices]
g_mag=mag[g_indices]
b_ra=ra[b_indices]
g_ra=ra[g_indices]


# In[15]:


plt.hist(b_mag, bins=100,color = bronze ) 
plt.hist(g_mag, bins=100,color = 'gold') 
plt.xlabel('mag')  
plt.ylabel('event')  
plt.title('distribution of mag Standard candle')
plt.xlim()
plt.show()


# In[18]:


b_coords = SkyCoord(ra=b_ra*u.degree, dec=b_dec*u.degree, frame='icrs')
g_coords = SkyCoord(ra=g_ra*u.degree, dec=g_dec*u.degree, frame='icrs')
# Plot the RA and Dec coordinates on a map
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='mollweide')
ax.scatter(b_coords.ra.wrap_at(180*u.degree).radian, b_coords.dec.radian, s=50, color=bronze)
ax.scatter(g_coords.ra.wrap_at(180*u.degree).radian, g_coords.dec.radian, s=50, color='gold')
# Add magnitude estimates as numbers for each point on the map
for i, (ra_deg, dec_deg, mag) in enumerate(zip(b_ra, b_dec, b_mag)):
    ax.text(np.radians(ra_deg), np.radians(dec_deg), f'{mag:.2f}',
            fontsize=8, color='black', ha='center', va='center')
for i, (ra_deg, dec_deg, mag) in enumerate(zip(g_ra, g_dec, g_mag)):
    ax.text(np.radians(ra_deg), np.radians(dec_deg), f'{mag:.2f}',
            fontsize=8, color='black', ha='center', va='center')
# Set plot title and labels
ax.set_title('Simulated events with magnitude on sky map 01')
ax.set_xlabel('RA (radian)')
ax.set_ylabel('Dec (radian)')

# Use astropy's MPL style
plt.style.use(astropy_mpl_style)


plt.show()


# In[24]:


headers = ["b_dec", "b_ra", "b_mag"]

# Combine the arrays into a list of lists
data = list(zip(b_dec, b_ra, b_mag))

# Print the table
print(tabulate(data, headers=headers))


# In[26]:


headers = ["g_dec", "g_ra", "g_mag"]

# Combine the arrays into a list of lists
data = list(zip(g_dec, g_ra, g_mag))

# Print the table
print(tabulate(data, headers=headers))


# In[ ]:




