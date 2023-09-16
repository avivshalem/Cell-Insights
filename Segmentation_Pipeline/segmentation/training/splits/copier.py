from pathlib import Path

import numpy as np

from random import shuffle
import shutil


path_train_sets = Path(r'C:\Users\avivs\PycharmProjects\env\TAU\Research\New-Segmentation\kit-sch-ge_2021_segmentation\training_sets')
cell_type = 'PhC-C2DL-PSC'
mode = 'GT'

img_ids = sorted((path_train_sets / "{}_{}".format(cell_type, mode)).glob('img*.png'))

img_ids_stem = []
for idx in img_ids:
    img_ids_stem.append(idx.stem.split('img_')[-1])

# Random 80%/20% split
shuffle(img_ids_stem)
train_ids = img_ids_stem[0:int(np.floor(0.8 * len(img_ids_stem)))]
val_ids = img_ids_stem[int(np.floor(0.8 * len(img_ids_stem))):]

train_val_ids = {'train': train_ids, 'val': val_ids}

for train_mode in ['train', 'val']:
    for idx in train_val_ids[train_mode]:
        source_path = path_train_sets / "{}_{}".format(cell_type, mode)
        target_path = path_train_sets / "{}_{}".format(cell_type, mode) / train_mode
        # if (source_path / "A" / ("img_{}.tif".format(idx))).exists():
        #     source_path = source_path / "A"
        # else:
        #     source_path = source_path / "B"
        shutil.copyfile(str(source_path / "img_{}.png".format(idx)),
                        str(target_path / "img_{}.png".format(idx)))
        shutil.copyfile(str(source_path / "dist_cell_{}.png".format(idx)),
                        str(target_path / "dist_cell_{}.png".format(idx)))
        shutil.copyfile(str(source_path / "dist_neighbor_{}.png".format(idx)),
                        str(target_path / "dist_neighbor_{}.png".format(idx)))
        shutil.copyfile(str(source_path / "mask_{}.png".format(idx)),
                        str(target_path / "mask_{}.png".format(idx)))