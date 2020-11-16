#import packages and functions
from pathlib import Path
from plot import plotFeatures
from arch import dame_arch_paths, dame_arch_data
#####################################################################################
# MANUAL INPUTS REQUIRED FROM USER

# change the path
path_to_folders = str(Path(r'C:\Users\Lisa Paillard\Desktop\Pour TFE\jHamON_data'))

# participants - adapt the list
part = ['01_Kevin', '02_rafaelopes', '03_charlesbarrand', '04_guilhem',\
        '05_leandre', '06_thomasmartine', '10_victor',\
        '11_youssouf', '12_sufyan', '16_julien', '34_nicolas']

# colors list for visualization: should be at least as long as part list
colors = [(1,0.26,0.11),(0.45,0.45,1),(0.8,0.2,1),(0,0.4,1),\
          (0,1,0.6),(0.2,0.8,0.2),(1,0.8,0),(0.8,0,0),\
          (1,0,0.4),(0.7,0.24,0),(0,0,0.8)]

#####################################################################################
# launch analysis automatically
archpaths = dame_arch_paths(path_to_folders=path_to_folders, participants = part)
arch_dict = dame_arch_data(archpaths)

# save results in dict object
from dictmanager import save_obj
save_obj(arch_dict, path_to_folders + '\\results')

# visualize results
plotFeatures(path_to_folders, '\\results', colors, part)