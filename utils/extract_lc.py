import urllib.request
import pandas as pd
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt

def read_from_file(file_path):
    lc = pd.read_csv(file_path)
    return lc

def read_from_open_catalog(name_='sn1987a',
                           allowed_bands= ['U', 'u', 'B', 'V', 'g', 'R', 'r', \
                                           'I', 'i', 'J', 'H', "u'", "g'", "r'", "i'"]
                                           ):

    if ' ' in name_:
        name_url = '%20'.join(name_.split(' '))
    else:
        name_url = name_
    url = 'https://api.astrocats.space/'+str(name_url)+'/photometry/time+magnitude+e_magnitude+band+instrument+u_time+source+upperlimit?format=csv'

    print ('Reading in '+ str(name_)+' ...')

    response = urllib.request.urlopen(url)
    cr = pd.read_csv(response)
    # bands = cr.groupby('band').size()

    # Get rid of data points with unknown bands
    print('Number of points removed because of Nan bands: ', np.sum(cr.band.isnull()))
    cr = cr[(cr.band.isnull()) != True]

    #Select data points that have a band that exists in bands
    print('Number of points removed because of undesired bands: ', np.sum(~cr.band.isin(allowed_bands)))
    cr = cr[cr.band.isin(allowed_bands)]

    # Get rid of data points with upper limit photometry
    print('Number of points removed because of being upperlimit: ', np.sum(cr.upperlimit == 'T'))
    cr = cr[cr.upperlimit != 'T']

    cr['magnitude'] = cr['magnitude'].apply(pd.to_numeric)
    cr['e_magnitude'] = cr['e_magnitude'].apply(pd.to_numeric)
    cr['time'] = cr['time'].apply(pd.to_numeric)

    cr["date"] = Time(cr['time'] + 2400000.5, format='jd').to_datetime()

    return cr

def plot_all_band(lc, name_, ax):
    print(np.unique(lc['band']))
    for i in np.unique(lc['band']):
        ax.errorbar(lc[lc.band == i].date, 
                    lc[lc.band == i].magnitude,
                    yerr=lc[lc.band == i].e_magnitude,
                    fmt ='o', label=i)
        
    ax.tick_params(axis='x', labelrotation=90)
    x = 'Plot of '+ str(name_) + ' in band '+str(i)
    plt.title(x)
    ax.set_xlabel('Time (MJD)')
    ax.set_ylabel('Magnitude')
    ax.invert_yaxis()
    plt.legend()

    