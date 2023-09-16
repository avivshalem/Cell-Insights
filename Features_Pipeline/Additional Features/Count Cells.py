import os
import tifffile as tiff
import cv2
import numpy as np
import pandas as pd
from skimage import measure

filename = 'output_number_cells'
path = r'\HCR'
path_list = os.listdir(path)
# dict_it = dict.fromkeys(['Experiment', dict.fromkeys(['Circularity', 'Circularity minimum enclosing circle', 'Circularity minimum enclosing ellipse',
#                          'Aspect Ratio', 'Eccentricity', 'Form Factor'])])
dict_it = {}

for folder in path_list:
    current_folder = os.path.join(path, folder)
    if not os.path.isdir(current_folder):
        continue
    dict_it[folder] = {'Number of Cells': dict.fromkeys(['start', 'end'])}

    current_folder = os.path.join(path, folder)
    folder_list = os.listdir(current_folder)

    startI = folder_list[0]
    endI = folder_list[-1]
    if len(folder_list) >= 164:
        endI = folder_list[163]

    start_image = os.path.join(current_folder, startI)
    read_image = tiff.imread(start_image)
    read_image_normed = np.round((read_image / np.max(read_image)) * 255).astype(np.uint8)
    labeled_image = measure.label(read_image_normed)
    start_n_objects = len(np.unique(labeled_image)) - 1

    end_image = os.path.join(current_folder, endI)
    read_image = tiff.imread(end_image)
    read_image_normed = np.round((read_image / np.max(read_image)) * 255).astype(np.uint8)
    labeled_image = measure.label(read_image_normed)
    end_n_objects = len(np.unique(labeled_image)) - 1

    for k, key in enumerate(dict_it[folder].keys()):
        dict_it[folder][key]['start'] = start_n_objects
        dict_it[folder][key]['end'] = end_n_objects
df = pd.DataFrame.from_dict(dict_it, orient='index')
df.to_excel(path + f'\{filename}.xlsx')

