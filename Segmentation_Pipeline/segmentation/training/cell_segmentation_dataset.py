import numpy as np
import tifffile as tiff
import imageio as imageio
from PIL import Image

from torch.utils.data import Dataset


class CellSegDataset(Dataset):
    """ Pytorch data set for instance cell nuclei segmentation """

    def __init__(self, root_dir, cell_type, gt_mode, mode='train', transform=lambda x: x):
        """

        :param root_dir: Directory containing all created training/validation data sets.
            :type root_dir: pathlib Path object.
        :param cell_type: cell type.
            :type cell_type: str
        :param gt_mode: Primary Track mode ('GT', 'ST', 'GT+ST', 'allGT', 'allST', 'allGT+allST')
            :type gt_mode: str
        :param mode: 'train' or 'val'.
            :type mode: str
        :param transform: transformss.
            :type transform:
        :return: Dict (image, cell_label, border_label, id).
        """
        # print(root_dir, cell_type, gt_mode)
        if "all" in gt_mode:
            self.img_ids = sorted((root_dir / gt_mode / mode).glob('img*.tif'))
        else:
            self.img_ids = sorted((root_dir / "{}_{}".format(cell_type, gt_mode) / mode).glob('img*.png'))
        self.mode = mode
        self.root_dir = root_dir
        self.transform = transform

    def __len__(self):
        # print(self.root_dir)
        return len(self.img_ids)

    def __getitem__(self, idx):

        img_id = self.img_ids[idx]

        img = np.asarray(Image.open(str(img_id))).astype(np.uint16)
        temp_shape = img.shape
        img = np.reshape(img,(temp_shape[0],temp_shape[1],1))

        dist_label_id = img_id.parent / ('dist_cell' + img_id.name.split('img')[-1])
        dist_neighbor_label_id = img_id.parent / ('dist_neighbor' + img_id.name.split('img')[-1])

        dist_label = np.asarray(Image.open(str(dist_label_id))).astype(np.float32)
        temp_shape = dist_label.shape
        dist_label = np.reshape(dist_label,(temp_shape[0],temp_shape[1],1))
        dist_neighbor_label = np.asarray(Image.open(str(dist_neighbor_label_id))).astype(np.float32)
        temp_shape = dist_neighbor_label.shape
        dist_neighbor_label = np.reshape(dist_neighbor_label,(temp_shape[0],temp_shape[1],1))

        sample = {'image': img,
                  'cell_label': dist_label,
                  'border_label': dist_neighbor_label,
                  'id': img_id.stem}

        sample = self.transform(sample)

        return sample
