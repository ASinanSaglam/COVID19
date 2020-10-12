import sys
sys.path.append('../../python-loader')
from pyMCDS import pyMCDS
import pickle, os

# mcds = pyMCDS("output{:08d}.xml".format(0), "output/")
# import IPython;IPython.embed()
# sys.exit()

if os.path.exists("virion_arr.pickle"):
    with open("virion_arr.pickle", "rb") as f:
        virion_dat = pickle.load(f)
else:
    virion_dat = []
    for i in range(0,145):
        mcds = pyMCDS("output{:08d}.xml".format(i), "output/")
        total_volume = 2494
        virion_dens = mcds.data['continuum_variables']['virion']['data'].sum()
        virions = virion_dens*total_volume
        virion_dat.append(virions)
    
    with open("virion_arr.pickle", "wb") as f:
        pickle.dump(virion_dat, f)

import matplotlib.pyplot as plt
import seaborn

plt.plot(virion_dat)
plt.xlabel("time (hours)")
plt.ylabel("total virion density * total volume")
plt.savefig("virions_2843cells_custom.png")
plt.close()
