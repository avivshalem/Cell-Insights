from pathlib import Path
import os

import numpy as np
from tifffile import imread

from tracker.export import ExportResults
from tracker.extract_data import get_img_files
from tracker.extract_data import get_indices_pandas
from tracker.tracking import TrackingConfig, MultiCellTracker
import pickle


class TrackArgs():
    def __init__(self, image_path, segmentation_path, results_path, delta_t, default_roi_size):
        self.img_path = image_path
        self.segm_path = segmentation_path
        self.res_path = results_path
        self.delta_t = delta_t
        self.default_roi_size = default_roi_size


def main(img_path, segm_path, res_path, delta_t=3, default_roi_size=2):

    args = TrackArgs(img_path, segm_path, res_path, delta_t, default_roi_size)

    img_path = Path(args.img_path)
    segm_path = Path(args.segm_path)
    res_path = Path(args.res_path)
    img_files = get_img_files(img_path)
    segm_files = get_img_files(segm_path, 'mask')

    # set roi size
    # assume img shape z,x,y
    dummy = np.squeeze(imread(segm_files[max(segm_files.keys())]))
    img_shape = dummy.shape
    masks = get_indices_pandas(imread(segm_files[max(segm_files.keys())]))
    m_shape = np.stack(masks.apply(lambda x: np.max(np.array(x), axis=-1) - np.min(np.array(x), axis=-1) +1))

    if len(img_shape) == 2:
        if len(masks) > 10:
            m_size = np.median(np.stack(m_shape)).astype(int)

            roi_size = tuple([m_size*args.default_roi_size, m_size*args.default_roi_size])
        else:
            roi_size = tuple((np.array(dummy.shape) // 10).astype(int))
    else:
        roi_size = tuple((np.median(np.stack(m_shape), axis=0) * args.default_roi_size).astype(int))

    config = TrackingConfig(img_files, segm_files, roi_size, delta_t=args.delta_t, cut_off_distance=None, allow_cell_division=False)
    tracker = MultiCellTracker(config)
    tracks = tracker()

    # Aviv - 26/7
    with open(os.path.join(args.res_path, 'tracks.pickle'), 'wb') as handle:
        pickle.dump(tracks, handle, protocol=pickle.HIGHEST_PROTOCOL)

    exporter = ExportResults()
    exporter(tracks, res_path, tracker.img_shape, time_steps=sorted(img_files.keys()))


if __name__ == '__main__':

    main()
