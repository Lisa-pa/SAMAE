from pathlib import Path

from arch import *

# manually change the path or ask for a user folder selection input
# path_to_folders = str(Path(r'C:\Users\Antonio\Dropbox\Paillard_Lisa\jHamON_data'))
path_to_folders = str(Path(r'C:\Users\Lisa Paillard\Desktop\Pour TFE\jHamON_data'))


archpaths = dame_arch_paths(path_to_folders=path_to_folders)
arch_dict = dame_arch_data(archpaths)


# # example of nested dicts 
# arch_dict['01_Kevin']['fam_1']['BF']['simple']['img_1']['architecture']['fsc_1'].keys()


#if path leads to panoramic image -> execute main_panoramic
#elif path leads to simple image -> execute main_simple