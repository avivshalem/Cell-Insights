# KIT-Sch-GE 2021 Segmentation - Simple Modify for Cell Insights

Segmentation method used for our submission to the 6th edition of the [ISBI Cell Tracking Challenge](http://celltrackingchallenge.net/) 2021 (Team KIT-Sch-GE).

## Prerequisites
* [Anaconda Distribution](https://www.anaconda.com/products/individual)
* A CUDA capable GPU
* Minimum / recommended RAM: 16 GiB / 32 GiB
* Minimum / recommended VRAM: 12 GiB / 24 GiB


### Data
```
python CreateDatasetBeforeTraining.py --train --cell_type "all" --mode "gt"
```
### Training
```
python CreateDatasetBeforeTraining.py --train --cell_type "all" --mode "gt"
```

## Publication ##
T. Scherr, K. Löffler, M. Böhland, and R. Mikut (2020). Cell Segmentation and Tracking using CNN-Based Distance Predictions and a Graph-Based Matching Strategy. PLoS ONE 15(12). DOI: [10.1371/journal.pone.0243219](https://doi.org/10.1371/journal.pone.0243219).

## License ##
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.