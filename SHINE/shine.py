#!/usr/bin/env python
# coding: utf-8

# Viewing solar data with SunPy
# ---
# In the first half of this tutorial, we'll look at how to search for, download, and plot remote sensing data using SunPy.

# Setting search parameters
# ---
# In order to search for some data, we have to select a time range and an instrument to search for. In addition, here we also specify the wavelength of interest.

# In[1]:


from sunpy.net import Fido, attrs
import astropy.units as u

timerange = attrs.Time('2018/10/31 16:00', '2018/10/31 17:00')
instrument = attrs.Instrument('AIA')
wavelength = attrs.Wavelength(193 * u.Angstrom)


# Searching for data
# ---
# Using the above defined search parameters, Fido can be used to search for data. For more information on searching for and downloading data see https://docs.sunpy.org/en/stable/guide/acquiring_data/index.html

# In[2]:


result = Fido.search(timerange, instrument, wavelength)
result


# Downloading data
# ---
# The results from a search can also be downloaded. The above 'results' object can be indexed like an array. The first index is for the data provider, and the second index selects the row of interest. In this example we just download the first iamge from the first provider, so use ``results[0, 0]``.
# 
# ``Fido.fetch`` returns a list of the local location of the downloaded files.

# In[3]:


downloaded_files = Fido.fetch(result[0, 0])
print(downloaded_files[0])


# Loading and plotting data
# ---
# Now we have downloaded some data, we can load it and plot it. ``sunpy.map.Map`` can be used to load any ``.fits`` file, creating a ``Map`` object. We can then take a look at the image stored by calling ``map.peek()``.

# In[4]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import sunpy.map

aiamap = sunpy.map.Map(downloaded_files[0])
ax = plt.subplot(projection=aiamap)
aiamap.plot(ax)


# Identifying a coronal hole
# ---
# In the above AIA 193 image, there is a dark patch in the middle of the disc. This is a coronal hole, which is the source of the fastest solar wind.
# 
# This observation can be connected to in-situ measurements of the solar wind at 1 AU by looking for the fast solar wind stream that emmenates from this coronal hole. Because the solar wind takes a finite amount of time to propagate from the Sun to Earth, we first do an order of magnitude estimate of this delay.

# Using astropy units
# ---
# To calculate the propagation delay we can use the ``astropy.units`` module. This provides an extension of normal numbers and arrays, and allows units to be attached. All the unit mathematics is calculated automatically, avoiding the need to keep track of specific units.

# In[5]:


import astropy.constants as const
# Take a typical fast solar wind speed of 500 km/s
vsw = 500 * u.km / u.s
# Assume the solar wind is relesed from the surface of the Sun, so
# the propagation distance is 1 AU - r_sun
d = (1 * u.au - const.R_sun)
# Calculate tne propagation time, and convert it to units of days
t = (d / vsw).to(u.d)
print(t)


# Downloading and importing in-situ data
# ---
# The ``heliopy.data`` module can be used to download and import a wide range of in situ datasets from various heliospheric missions. In this example we use data from OMNI, which provides measurements of the solar wind at the orbit of the Earth.

# In[6]:


from heliopy.data import omni
from datetime import datetime


# In[7]:


starttime = datetime(2018, 10, 1)
endtime = datetime(2018, 12, 1)
data = omni.low(starttime, endtime)


# The data is stored in the ``data`` object. We can print the available columns in this object:

# In[8]:


for col in data.columns:
    print(col)


# Plotting in-situ data
# ---
# Matplotlib can be used to plot the downloaded data. In this example we plot the solar wind speed and the magnetic field clock angle, to see different polarity solar wind streams.
# 
# We also add a vertical line where the stream is fast stream is predicted to have arrived using the above back-of-the-envelope calculate. We can see that it lines up nicely with a fast solar wind stream that has speeds of 500 - 600 km/s.

# In[9]:


from astropy.visualization import quantity_support
quantity_support()

fig, axs = plt.subplots(figsize=(10, 6), nrows=3, sharex=True)

ax = axs[0]
ax.plot(data.index, data.quantity('Plasma Flow Speed'), label='$v_{sw}$')
ax.axvline(datetime(2018, 11, 5), color='k')

ax = axs[1]
ax.plot(data.index, data.quantity('|B|'), label='$|B|$')

ax = axs[2]
ax.plot(data.index, data.quantity('Proton Density'), label='$n_{p}$')


# Improving figure formatting
# ---

# In[10]:


# Add a legend to each axes
for ax in axs:
    ax.legend()
# Make the x-axis formatting nicer
fig.autofmt_xdate()
fig.subplots_adjust(hspace=0)
# Move the middle y-axis to the right hand side
axs[1].yaxis.tick_right()
# Line up the y-axis labels
fig.align_labels()

# Show figure again
fig


# Save a copy of the figure
# ---

# In[11]:


fig.savefig('tseries.pdf', bbox_inches='tight')

