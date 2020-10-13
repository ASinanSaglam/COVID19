import sys
sys.path.append('../../python-loader')
from pyMCDS import pyMCDS
import pickle, os

# mcds = pyMCDS("output{:08d}.xml".format(0), "output/")
# import IPython;IPython.embed()
# sys.exit()

if os.path.exists("data_arr.pickle"):
    with open("data_arr.pickle", "rb") as f:
        dat = pickle.load(f)
    virion_dat, disc_vir, cell_dat, immune_dat, cont_immune = dat
else:
    virion_dat = []
    cell_dat = []
    disc_vir = []
    immune_dat = []
    cont_immune = []
    for i in range(0,101):
        mcds = pyMCDS("output{:08d}.xml".format(i), "output/")
        # import IPython;IPython.embed()
        # sys.exit() 
        total_volume = 2494
        virion_dens = mcds.data['continuum_variables']['virion']['data'].sum()
        virions = virion_dens*total_volume
        virion_dat.append(virions)

        disc_vir_cnt = mcds.data['discrete_cells']['virion'].sum()
        disc_vir.append(disc_vir_cnt)

        live_cell_cnt = (mcds.data['discrete_cells']['cell_type']==1).sum()
        cell_dat.append(live_cell_cnt)

        immune_cnt = (mcds.data['discrete_cells']['cell_type']>1).sum()
        immune_dat.append(immune_cnt)

        cont_immune_cnt = mcds.data['continuum_variables']['interferon 1']['data'].sum()
        cont_immune_cnt += mcds.data['continuum_variables']['pro-inflammatory cytokine']['data'].sum()
        cont_immune_cnt += mcds.data['continuum_variables']['chemokine']['data'].sum()
        cont_immune_cnt += mcds.data['continuum_variables']['debris']['data'].sum()
        cont_immune_cnt = cont_immune_cnt * total_volume
        cont_immune.append(cont_immune_cnt)
        
    dat = (virion_dat, disc_vir, cell_dat, immune_dat, cont_immune)
    
    with open("data_arr.pickle", "wb") as f:
        pickle.dump(dat, f)

import matplotlib.pyplot as plt
import seaborn

plt.plot(virion_dat, label="continuum virions")
plt.plot(disc_vir, label="discrete virions")
plt.plot(cell_dat, label="live cells")
plt.plot(immune_dat, label="immune cells")
plt.plot(cont_immune, label="continuum immune resp")
plt.legend(frameon=False)
plt.gca().set(yscale="log")
plt.xlabel("time (hours)")
plt.ylabel("count (or total virion density * total volume)")
plt.savefig("custom_viral_inf.png")
plt.close()
