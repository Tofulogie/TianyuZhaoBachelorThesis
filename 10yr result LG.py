#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization import astropy_mpl_style


# In[3]:


gfile = np.genfromtxt('10yr gold LG.txt')
bfile = np.genfromtxt('10yr bronze LG.txt')

g_mag=gfile[:, 2]
b_mag=bfile[:, 2]
g_ra=gfile[:, 1]
b_ra=bfile[:, 1]
g_dec=gfile[:, 0]
b_dec=bfile[:, 0]

x_value = 19.5 
plt.axvline(x=x_value, color='r', linestyle='--')
x1_value = 22  
plt.axvline(x=x1_value, color='r', linestyle='--')
plt.hist(b_mag, bins=100) 
plt.xlabel('magnitude')  
plt.ylabel('Event') 
plt.xlim()
plt.show()


# In[8]:


x_value = 19.5  
plt.axvline(x=x_value, color='r', linestyle='--')
x1_value = 22  
plt.axvline(x=x1_value, color='r', linestyle='--')
plt.hist(g_mag, bins=100) 
plt.xlabel('magnitude')  
plt.ylabel('Event')   
plt.xlim()
plt.show()


# In[9]:


plt.hist(g_dec, bins=100) 
plt.xlabel('dec')  
plt.ylabel('Event')  
plt.xlim(-90,90)
plt.show()


# In[10]:


plt.hist(b_dec, bins=100) 
plt.xlabel('dec')  
plt.ylabel('Event')  
plt.xlim(-90,90)
plt.show()


# In[18]:


bronze = (0.8, 0.5, 0.2) 
b_coords = SkyCoord(ra=b_ra*u.degree, dec=b_dec*u.degree, frame='icrs')
g_coords = SkyCoord(ra=g_ra*u.degree, dec=g_dec*u.degree, frame='icrs')
# Plot the RA and Dec coordinates on a map
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='mollweide')
ax.scatter(b_coords.ra.wrap_at(180*u.degree).radian, b_coords.dec.radian, s=50, color='blue',alpha=0.5)
ax.scatter(g_coords.ra.wrap_at(180*u.degree).radian, g_coords.dec.radian, s=50, color='red',alpha=0.5)
# Add magnitude estimates as numbers for each point on the map
#for i, (ra_deg, dec_deg, mag) in enumerate(zip(b_ra, b_dec, b_mag)):
#    ax.text(np.radians(ra_deg), np.radians(dec_deg), f'{mag:.2f}',
#            fontsize=5, color='black', ha='center', va='center')
#for i, (ra_deg, dec_deg, mag) in enumerate(zip(g_ra, g_dec, g_mag)):
#   ax.text(np.radians(ra_deg), np.radians(dec_deg), f'{mag:.2f}',
#            fontsize=5, color='black', ha='center', va='center')
# Set plot title and labels

ax.set_xlabel('RA (radian)')
ax.set_ylabel('Dec (radian)')

# Use astropy's MPL style
plt.style.use(astropy_mpl_style)


plt.show()


# In[12]:


condition_lower = 19.5 < b_mag
condition_upper = b_mag <= 22

# Combine conditions using the bitwise AND operator
combined_condition = condition_lower & condition_upper

b0_indices = np.where(b_mag<= 19.5)
# Use np.where to get the indices where the combined condition is True
b1_indices = np.where(combined_condition)

condition_lower_g = 19.5 < g_mag
condition_upper_g = g_mag <= 22

# Combine conditions using the bitwise AND operator
combined_condition_g = condition_lower_g & condition_upper_g

g0_indices = np.where(g_mag<= 19.5)
# Use np.where to get the indices where the combined condition is True
g1_indices = np.where(combined_condition_g)


# In[19]:


bronze = (0.8, 0.5, 0.2) 
b_coords = SkyCoord(ra=b_ra[b0_indices]*u.degree, dec=b_dec[b0_indices]*u.degree, frame='icrs')
g_coords = SkyCoord(ra=g_ra[g0_indices]*u.degree, dec=g_dec[g0_indices]*u.degree, frame='icrs')
# Plot the RA and Dec coordinates on a map
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='mollweide')
ax.scatter(b_coords.ra.wrap_at(180*u.degree).radian, b_coords.dec.radian, s=50, color='blue',alpha=0.5)
ax.scatter(g_coords.ra.wrap_at(180*u.degree).radian, g_coords.dec.radian, s=100, color='red',alpha=0.5)
# Add magnitude estimates as numbers for each point on the map
for i, (ra_deg, dec_deg, mag) in enumerate(zip(b_ra[b0_indices], b_dec[b0_indices], b_mag[b0_indices])):
    ax.text(np.radians(ra_deg), np.radians(dec_deg), f'{mag:.2f}',
            fontsize=8, color='black', ha='center', va='center')
for i, (ra_deg, dec_deg, mag) in enumerate(zip(g_ra[g0_indices], g_dec[g0_indices], g_mag[g0_indices])):
    ax.text(np.radians(ra_deg), np.radians(dec_deg), f'{mag:.2f}',
            fontsize=8, color='black', ha='center', va='center')
# Set plot title and labels

ax.set_xlabel('RA (radian)')
ax.set_ylabel('Dec (radian)')

# Use astropy's MPL style
plt.style.use(astropy_mpl_style)


plt.show()


# In[20]:


bronze = (0.8, 0.5, 0.2) 
b_coords = SkyCoord(ra=b_ra[b1_indices]*u.degree, dec=b_dec[b1_indices]*u.degree, frame='icrs')
g_coords = SkyCoord(ra=g_ra[g1_indices]*u.degree, dec=g_dec[g1_indices]*u.degree, frame='icrs')
# Plot the RA and Dec coordinates on a map
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='mollweide')
ax.scatter(b_coords.ra.wrap_at(180*u.degree).radian, b_coords.dec.radian, s=50, color='blue',alpha=0.5)
ax.scatter(g_coords.ra.wrap_at(180*u.degree).radian, g_coords.dec.radian, s=100, color='red',alpha=0.5)
# Add magnitude estimates as numbers for each point on the map
for i, (ra_deg, dec_deg, mag) in enumerate(zip(b_ra[b1_indices], b_dec[b1_indices], b_mag[b1_indices])):
    ax.text(np.radians(ra_deg), np.radians(dec_deg), f'{mag:.2f}',
            fontsize=8, color='black', ha='center', va='center')
for i, (ra_deg, dec_deg, mag) in enumerate(zip(g_ra[g1_indices], g_dec[g1_indices], g_mag[g1_indices])):
    ax.text(np.radians(ra_deg), np.radians(dec_deg), f'{mag:.2f}',
           fontsize=8, color='black', ha='center', va='center')
# Set plot title and labels

ax.set_xlabel('RA (radian)')
ax.set_ylabel('Dec (radian)')

# Use astropy's MPL style
plt.style.use(astropy_mpl_style)


plt.show()


# In[ ]:





# In[ ]:




