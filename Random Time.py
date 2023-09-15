#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from astropy.time import Time, TimeDelta
import numpy as np

# Define a start and end date range
start_date = Time("2024-01-01")
end_date = Time("2024-12-31")

# Generate Time within the date range
random_date_offsets = TimeDelta(np.random.uniform(0, 1, 2) * (end_date - start_date), format="jd")
random_times = start_date + random_date_offsets

print(random_times)

