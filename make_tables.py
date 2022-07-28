
import numpy as np


import pickle

from astropy.io import ascii
from astropy.table import vstack, Table

import flares_utility.table as tab



import flares_utility.analyse

# ----------------------------------------------------------------------
# --- open data and load analyser

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed_v2.hdf5'
a = flares_utility.analyse.analyse_flares(filename, default_tags = False)


# tags = ['005_z010p000','006_z009p000','007_z008p000','008_z007p000','009_z006p000','010_z005p000']
# zeds = np.array([float(tag[5:].replace('p','.')) for tag in tags])

tags = ['000_z015p000','001_z014p000','002_z013p000','003_z012p000','004_z011p000','005_z010p000','006_z009p000','007_z008p000','008_z007p000','009_z006p000','010_z005p000']
zeds = np.array([float(tag[5:].replace('p','.')) for tag in tags])

print(zeds)




# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]


ages = pickle.load(open('../analysis/sf/moments_and_percentiles.p','rb'))
# metallicities = pickle.load(open('data/HOFe.p','rb'))

quantities = []

quantities.append({'path': 'Galaxy/Metallicity/', 'dataset': f'MassWeightedStellarZ', 'name': 'Zstar', 'log10': True})
quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar', 'log10': True})
quantities.append({'path': f'Galaxy/SFR_aperture/30', 'dataset': f'50Myr', 'name': f'SFR', 'log10': True})


x = 'log10Mstar'

if x[:5] == 'log10':
    x_ = x[5:]
else:
    x_ = x


for y in ['log10Zstar', 'log10SFR', 'log10sSFR', 'age']:

    if y[:5] == 'log10':
        y_ = y[5:]
    else:
        y_ = y


    D = {}
    s = {}
    tables = []

    for tag, z in zip(tags, zeds):

        # --- get quantities (and weights and deltas)
        D = a.get_datasets(tag, quantities)
        D['log10sSFR'] = D['log10SFR'] - D['log10Mstar'] + 9
        D['age'] = ages[z]['P0.5']

        t = tab.binned_single_z(D, x, y, bins = np.arange(8., 12.0, 0.2))

        t['z'] = z*np.ones(len(t[x].data))

        tables.append(t)


    table = vstack(tables)

    table.meta['name'] = 'FLARES'
    table.meta['redshifts'] = list(set(table['z'].data))
    table.meta['references'] = ['2021MNRAS.500.2127L']
    table.meta['x'] = x
    table.meta['y'] = y

    table.write(f'data/{x_}/{y_}/flares.ecsv', overwrite = True)
    table.write(f'/Users/stephenwilkins/Dropbox/Research/modules/flare_data/flags_data/data/ScalingRelations/{x_}/{y_}/models/flares.ecsv', overwrite = True)
