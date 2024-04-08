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

        # Make directory for the trained models
        path_models.mkdir(exist_ok=True)

        # print('Create training sets for all cell types ...')
        # create_ctc_training_sets(path_data=path_datasets,
        #                          path_train_sets=path_train_data,
        #                          cell_types=paths['cell_types'])

        for cell_type in cell_types:

            for architecture in settings['methods']:

                # Get model names and how many iterations/models need to be trained
                if "all" in args.mode:
                    model_name = '{}_{}'.format(args.mode, architecture[2])
                else:
                    model_name = '{}_{}_{}'.format(cell_type, args.mode, architecture[2])
                if args.mode == 'GT' and architecture[-1]:  # auto-encoder pre-training only for GT
                    model_name += '-auto'
                    if cell_type in ['Fluo-N3DH-SIM+', 'Fluo-N2DH-SIM+']:  # not needed for simulated data
                        continue
                num_trained_models = len(list(path_models.glob('{}_model*.pth'.format(model_name))))
                if "all" in args.mode or "ST" in args.mode:
                    iterations = settings['iterations'] - num_trained_models
                else:
                    iterations = settings['iterations_GT_single_celltype'] - num_trained_models

                if iterations <= 0:
                    continue

                for i in range(iterations):  # Train multiple models

                    run_name = utils.unique_path(path_models, model_name + '_model_{:02d}.pth').stem

                    train_configs = {'architecture': architecture[0],
                                     'batch_size': settings['batch_size'],
                                     'batch_size_auto': settings['batch_size_auto'],
                                     'label_type': architecture[1],
                                     'loss': architecture[3],
                                     'num_gpus': num_gpus,
                                     'optimizer': architecture[2],
                                     'run_name': run_name
                                     }

                    net = unets.build_unet(unet_type=train_configs['architecture'][0],
                                           act_fun=train_configs['architecture'][2],
                                           pool_method=train_configs['architecture'][1],
                                           normalization=train_configs['architecture'][3],
                                           device=device,
                                           num_gpus=num_gpus,
                                           ch_in=1,
                                           ch_out=1,
                                           filters=train_configs['architecture'][4])

                    if 'auto' in run_name:  # Pre-training of the Encoder
                        net_auto = unets.build_unet(unet_type='AutoU',
                                                    act_fun=train_configs['architecture'][2],
                                                    pool_method=train_configs['architecture'][1],
                                                    normalization=train_configs['architecture'][3],
                                                    device=device,
                                                    num_gpus=num_gpus,
                                                    ch_in=1,
                                                    ch_out=1,
                                                    filters=train_configs['architecture'][4])

                        data_transforms_auto = augmentors(label_type='auto', min_value=0, max_value=65535)

                        # Load training and validation set
                        datasets = AutoEncoderDataset(root_dir=path_datasets,
                                                      cell_type=cell_type,
                                                      gt_train_dir=path_results / 'training_sets',
                                                      transform=data_transforms_auto)

                        # Train model
                        train_auto(net=net_auto,
                                   dataset=datasets,
                                   configs=train_configs,
                                   device=device,
                                   path_models=path_models)

                        # Load weights
                        if num_gpus > 1:
                            net_auto.module.load_state_dict(torch.load(str(path_models / '{}.pth'.format(run_name)),
                                                                       map_location=device))
                        else:
                            net_auto.load_state_dict(torch.load(str(path_models / '{}.pth'.format(run_name)),
                                                                map_location=device))

                        pretrained_dict = net_auto.state_dict()
                        net_dict = net.state_dict()

                        # 1. filter out unnecessary keys
                        pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in net_dict}
                        # 2. overwrite entries in the existing state dict
                        net_dict.update(pretrained_dict)
                        # 3. load the new state dict
                        net.load_state_dict(net_dict)

                        del net_auto

                    # The training images are uint16 crops of a min-max normalized image
                    data_transforms = augmentors(label_type=train_configs['label_type'], min_value=0, max_value=65535)
                    train_configs['data_transforms'] = str(data_transforms)

                    # Load training and validation set
                    print(path_results, cell_type, args.mode)
                    datasets = {x: CellSegDataset(root_dir=path_results / 'training_sets',
                                                  mode=x,
                                                  cell_type=cell_type,
                                                  gt_mode=args.mode,
                                                  transform=data_transforms[x])
                                for x in ['train', 'val']}

                    # Get number of training epochs depending on dataset size (just roughly to decrease training time):
                    if len(datasets['train']) + len(datasets['val']) >= 1000:
                        print('>1000')
                        train_configs['max_epochs'] = 200
                    elif len(datasets['train']) + len(datasets['val']) >= 500:
                        train_configs['max_epochs'] = 240
                    elif len(datasets['train']) + len(datasets['val']) >= 200:
                        train_configs['max_epochs'] = 320
                    elif len(datasets['train']) + len(datasets['val']) >= 100:
                        train_configs['max_epochs'] = 400
                    elif len(datasets['train']) + len(datasets['val']) >= 50:
                        train_configs['max_epochs'] = 480
                    else:
                        train_configs['max_epochs'] = 560

                    # Write information to json-file
                    utils.write_train_info(configs=train_configs, path=path_models)

                    # Train model
                    best_loss = train(net=net,
                                      datasets=datasets,
                                      configs=train_configs,
                                      device=device,
                                      path_models=path_models)

                    if train_configs['optimizer'] == 'ranger':  # Fine-tune with cosine annealing
                        net = unets.build_unet(unet_type=train_configs['architecture'][0],
                                               act_fun=train_configs['architecture'][2],
                                               pool_method=train_configs['architecture'][1],
                                               normalization=train_configs['architecture'][3],
                                               device=device,
                                               num_gpus=num_gpus,
                                               ch_in=1,
                                               ch_out=1,
                                               filters=train_configs['architecture'][4])

                        # Get best weights as starting point
                        if num_gpus > 1:
                            net.module.load_state_dict(torch.load(str(path_models / '{}.pth'.format(run_name)),
                                                                  map_location=device))
                        else:
                            net.load_state_dict(torch.load(str(path_models / '{}.pth'.format(run_name)),
                                                           map_location=device))
                        _ = train(net=net,
                                  datasets=datasets,
                                  configs=train_configs,
                                  device=device,
                                  path_models=path_models,
                                  best_loss=best_loss)


if __name__ == "__main__":

    main()