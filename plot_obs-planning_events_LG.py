#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import astropy_mpl_style, quantity_support

plt.style.use(astropy_mpl_style)
quantity_support()
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time, TimeDelta
from astropy import coordinates as coord
from astropy.io import ascii


# In[2]:


file = np.genfromtxt('LG_Neutrinos+random Mag.txt')
typpe = file[:,0]
dec_deg = file[:,1]
ra_deg = file[:,2]
mag = file[:,3]
# Create a SkyCoord object
sky_coord = coord.SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)

with open('LG_Neutrinos+random Mag.txt', 'r') as file:
    lines = file.readlines()
    times = [line.split()[4] for line in lines]

time_objects = Time(times, format='isot', scale="utc")
time_objects = time_objects+ TimeDelta(np.random.uniform(0, 1), format="jd") * 1

LAST = EarthLocation(lat= 30.0529838*u.deg, lon= 35.0407331*u.deg, height=415*u.m)
#DOI 10.1088/1538-3873/acd8f0

utcoffset = 3*u.hour  
times = Time(time_objects) - utcoffset
sky_coordaltaz = sky_coord.transform_to(AltAz(obstime=times,location=LAST))


# In[3]:


for i in range(32):
    midnights = Time(times) - utcoffset
    delta_midnights = np.linspace(-2, 10, 100)*u.hour

    frame_1stnights = AltAz(obstime=midnights[i]+delta_midnights,location=LAST)
    sky_coordaltazs_1stnight = sky_coord[i].transform_to(frame_1stnights)
    sky_coordairmasss_1stnight = sky_coordaltazs_1stnight.secz

    from astropy.coordinates import get_sun

    delta_midnights = np.linspace(15*24, 25*24, 1000)*u.hour #15-25days after neutrinos events



    times_15days = midnights[i] + delta_midnights
    frame_15days = AltAz(obstime=times_15days, location=LAST)
    sunaltazs_15days = get_sun(times_15days).transform_to(frame_15days)

    from astropy.coordinates import get_body

    moon_15days = get_body("moon", times_15days)
    moonaltazs_15days = moon_15days.transform_to(frame_15days)

    sky_coordaltazs_15days = sky_coord[i].transform_to(frame_15days)

    plt.plot(delta_midnights, sunaltazs_15days.alt, color='r', label='Sun')
    plt.plot(delta_midnights, moonaltazs_15days.alt, color=[0.75]*3, ls='--', label='Moon')
    plt.scatter(delta_midnights, sky_coordaltazs_15days.alt,
                c=sky_coordaltazs_15days.az, label=f'Event{i+1}', lw=0, s=8,
                cmap='viridis')
    plt.fill_between(delta_midnights, 0*u.deg, 90*u.deg,
                     sunaltazs_15days.alt < -0*u.deg, color='0.5', zorder=0)
    plt.fill_between(delta_midnights, 0*u.deg, 90*u.deg,
                     sunaltazs_15days.alt < -18*u.deg, color='k', zorder=0)
    plt.colorbar().set_label('Azimuth [deg]')
    plt.legend(loc='upper left')
    plt.xlim(24*15*u.hour, 24*16*u.hour)
    #plt.xticks((np.arange(13)*2-12)*u.hour)
    plt.ylim(0*u.deg, 90*u.deg)
    plt.xlabel('The time after Neutrino Alerts[hour]')
    plt.ylabel('Altitude [deg]')
    plt.title('Visibility in Israel Midnight')
    plt.show()
    print(f'Event{i+1}',max(sky_coordaltazs_15days.alt[sunaltazs_15days.alt < -18*u.deg]))

