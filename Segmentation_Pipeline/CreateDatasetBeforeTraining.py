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

def main():

    random.seed()
    np.random.seed()

    # Get arguments
    parser = argparse.ArgumentParser(description='Cell Segmentation')
    parser.add_argument('--train', '-t', default=True, action='store_true', help='Train new models')
    parser.add_argument('--evaluate', '-e', default=False, action='store_true', help='Evaluate models')
    parser.add_argument('--inference', '-i', default=False, action='store_true', help='Inference')
    parser.add_argument('--save_raw_pred', '-s', default=False, action='store_true', help='Save raw predictions')
    parser.add_argument('--cell_type', '-c', default='all', type=str, help='Cell type')
    parser.add_argument('--mode', '-m', default='GT', type=str, help='Mode for training')
    parser.add_argument('--th_seed', '-th_s', default=0.45, type=float, help='Seed threshold')
    parser.add_argument('--th_cell', '-th_c', default=0.08, type=float, help='Cell size threshold')
    parser.add_argument('--apply_clahe', '-ac', default=False, action='store_true', help='Apply CLAHE')
    parser.add_argument('--n_splitting', '-ns', default=40, type=int, help='Cell number to apply local splitting (only 3D)')
    parser.add_argument('--scale', '-sc', default=1.0, type=float, help='Scale for down-/upsampling (inference)')
    parser.add_argument('--batch_size', '-bs', default=8, type=int, help='Batch size (inference)')
    parser.add_argument('--multi_gpu', '-mgpu', default=True, action='store_true', help='Use multiple GPUs')
    parser.add_argument('--artifact_correction', default=False, action='store_true', help='Artifact correction (only for very dense cells, e.g., HSC')
    parser.add_argument('--fuse_z_seeds', default=False, action='store_true', help='Fuse seeds in z-direction')
    parser.add_argument('--cuda', default=True, action='store_true', help='Use CUDA')

    args = parser.parse_args()

    with open(Path.cwd() / 'cell_segmentation_train_settings.json') as f:
        settings = json.load(f)

    paths = {
              "cell_types":
            [
                "PhC-C2DL-PSC"
            ],
              "path_ctc_metric": "EvaluationSoftware/",
              "path_data": "data/",
              "path_results": "results/"
            }

    # Paths
    path_datasets = Path(paths['path_data'])
    path_results = Path(paths['path_results'])
    if path_results == '':
        path_results = path_datasets
    path_models = path_results / 'segmentation_models'
    path_best_models = path_datasets / 'Models' / 'Model'
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
    print('using ', device)
    if str(device) == 'cuda':
        torch.backends.cudnn.benchmark = True
    if args.multi_gpu:
        num_gpus = torch.cuda.device_count()
    else:
        num_gpus = 1
    print('multi gpus ', num_gpus)

    if args.train:  # Train model from scratch

        print('Create training sets for all cell types ...')
        create_ctc_training_sets(path_data=path_datasets,
                                 path_train_sets=path_train_data,
                                 cell_types=paths['cell_types'])


if __name__ == "__main__":

    main()