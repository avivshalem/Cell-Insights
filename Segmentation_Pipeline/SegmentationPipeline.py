import argparse
import json
import numpy as np
import random
import torch
import warnings

from pathlib import Path

from segmentation.inference.inference import inference_2d_ctc, inference_3d_ctc
from segmentation.training.cell_segmentation_dataset import CellSegDataset
from segmentation.training.autoencoder_dataset import AutoEncoderDataset
from segmentation.training.create_training_sets import create_ctc_training_sets, create_sim_training_sets
from segmentation.training.mytransforms import augmentors
from segmentation.training.training import train, train_auto
from segmentation.utils import utils, unets
from segmentation.utils.metrics import ctc_metrics, count_det_errors

warnings.filterwarnings("ignore", category=UserWarning)


class SegmentArgs():
    def __init__(self, th_seed, th_cell, apply_clahe, savePath, resultPath, artifact_correction, batch_size, cuda, multi_gpu):
        self.th_seed = th_seed
        self.th_cell = th_cell
        self.apply_clahe = apply_clahe
        self.savePath = savePath
        self.resultPath = resultPath
        self.artifact_correction = artifact_correction
        self.batch_size = batch_size
        self.cuda = cuda
        self.multi_gpu = multi_gpu
        self.cell_type = 'all'
        self.mode = 'GT'
        self.scale = 1.0
        self.save_raw_pred = False
        self.n_splitting = 40
        self.fuse_z_seeds = False


def main(th_seed, th_cell, apply_clahe, savePath, resultPath, artifact_correction=False, batch_size=1, cuda=True, multi_gpu=False):

    random.seed()
    np.random.seed()

    # Get arguments
    args = SegmentArgs(th_seed, th_cell, apply_clahe, savePath, resultPath, artifact_correction, batch_size, cuda, multi_gpu)

    paths = {
              "cell_types":
            [
                "PhC-C2DL-PSC"
            ],
              "path_ctc_metric": "C:/Users/avivs/PycharmProjects/env/TAU/Research/New-Segmentation/EvaluationSoftware/",
              "path_data": "C:/Users/avivs/PycharmProjects/env/TAU/Research/New-Segmentation/cell_tracking_challenge/",
              "path_results": "C:/Users/avivs/PycharmProjects/env/TAU/Research/New-Segmentation/kit-sch-ge_2021_segmentation/"
            }

    # Paths
    path_datasets = Path(paths['path_data'])
    path_results = Path(paths['path_results'])
    if path_results == '':
        path_results = path_datasets
    path_models = path_results / 'segmentation_models'
    path_best_models = path_datasets / 'kit-sch-ge_2021' / 'SW'
    path_train_data = path_results / 'training_sets'
    path_ctc_metric = Path(paths['path_ctc_metric'])
    if args.cell_type == 'all':
        cell_types = paths['cell_types']
    else:
        cell_types = [args.cell_type]

    # Set device for using CPU or GPU
    if args.cuda:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = "cpu"
    # print('using ', device)
    if str(device) == 'cuda':
        torch.backends.cudnn.benchmark = True
        if args.multi_gpu:
            num_gpus = torch.cuda.device_count()
        else:
            num_gpus = 1
    elif str(device) == 'cpu':
        num_gpus = 1
    # print('multi gpus ', num_gpus)


    for cell_type in cell_types:

        if "all" in args.mode:
            model = (path_best_models / "{}_model.pth".format(args.mode))
        else:
            model = (path_best_models / "{}_{}_model.pth".format(cell_type, args.mode))


        path_seg_results = path_datasets / 'challenge_datasets' / cell_type / 'KIT-Sch-GE_2021' / args.mode\
                           / 'CSB' / "RES"
        path_seg_results.mkdir(parents=True, exist_ok=True)
        print('Inference using {} on {}: th_seed: {}, th_cell: {}, scale: {}, clahe:{}, cuda:{}'.format(model.stem,
                                                                                                  cell_type,
                                                                                                  args.th_seed,
                                                                                                  args.th_cell,
                                                                                                  args.scale,
                                                                                                  args.apply_clahe,
                                                                                                  str(device)))

        if args.fuse_z_seeds:
            print('Seed fusion in z-direction: on')


        if '2D' in cell_type:
            # print('2d')
            inference_2d_ctc(model=model,
                             data_path=Path(args.savePath),
                             result_path=Path(args.resultPath),
                             device=device,
                             batchsize=args.batch_size,
                             args=args,
                             num_gpus=num_gpus)
        else:
            inference_3d_ctc(model=model,
                             data_path=path_datasets / 'challenge_datasets' / cell_type / '01',
                             result_path=path_seg_results,
                             device=device,
                             batchsize=args.batch_size,
                             args=args,
                             num_gpus=num_gpus)


if __name__ == "__main__":

    main()