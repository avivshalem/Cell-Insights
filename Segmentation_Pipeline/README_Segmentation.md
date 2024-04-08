# KIT-Sch-GE 2021 Segmentation - Simple Modify for Cell Insights

Segmentation method used for our submission to the 6th edition of the [ISBI Cell Tracking Challenge](http://celltrackingchallenge.net/) 2021 (Team KIT-Sch-GE).

## Prerequisites
* [Anaconda Distribution](https://www.anaconda.com/products/individual)
* A CUDA capable GPU
* Minimum / recommended RAM: 16 GiB / 32 GiB
* Minimum / recommended VRAM: 12 GiB / 24 GiB


### Download Dataset
1. Download the datasets from: https://drive.google.com/file/d/1BqgVUxOMg7cNq_7u7oqaBP6ePGAQYoBI/view?usp=drive_link
2. Extract the two folders inside in this current directory

### Data Augmentation
```
Augmentations.py
```

### Data Preparation
1. Move the data as follows:
   2. Augmented data to: ..\data\training_datasets\PhC-C2DL-PSC\01
   3. Labeled augmented data to: ..\data\training_datasets\PhC-C2DL-PSC\01_GT\SEG
   4. Alternatively, use the example dataset. After extraction it is structured correctly
2. Run:
```
python CreateDatasetBeforeTraining.py --train --cell_type "all" --mode "gt"
```
3. The final dataset is created in the Results folder 
### Training
1. Run:
```
python Train.py --train --cell_type "all" --mode "gt"
```
2. The model are created in the Results folder
## Publication ##
T. Scherr, K. Löffler, M. Böhland, and R. Mikut (2020). Cell Segmentation and Tracking using CNN-Based Distance Predictions and a Graph-Based Matching Strategy. PLoS ONE 15(12). DOI: [10.1371/journal.pone.0243219](https://doi.org/10.1371/journal.pone.0243219).

## License ##
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.