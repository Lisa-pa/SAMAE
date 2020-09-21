from pathlib import Path

from arch import dame_arch_paths, dame_arch_data

# manually change the path or ask for a user folder selection input
# path_to_folders = str(Path(r'C:\Users\Antonio\Dropbox\Paillard_Lisa\jHamON_data'))
path_to_folders = str(Path(r'C:/Users/Lisa Paillard/Desktop/Pour TFE/jHamON_data2'))


archpaths = dame_arch_paths(path_to_folders=path_to_folders)
arch_dict = dame_arch_data(archpaths)

from dictmanager import save_obj
save_obj(arch_dict, path_to_folders + '\\results')

# # example of nested dicts 
# arch_dict['01_Kevin']['fam_1']['BF']['simple']['img_1']['architecture']['fsc_1'].keys()