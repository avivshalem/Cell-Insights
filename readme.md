# Cell-Insights

_Cell-Insights_ are the developed tools used in paper:

_Aberrant migration features in primary skin fibroblasts of Huntington’s disease patients holds potential for unraveling disease progression using an image based machine learning tool._
Saja Gharaba, Aviv Shalem, et. al.

## Main Usage - The Full Tracker Pipeline
This pipeline will run segmentation, tracking, and x-y coordinates extractions

1. Download the models and example data from: https://drive.google.com/file/d/1JphZWk8T7Hbekd_bB0vGgoV39_BE9tOb/view?usp=sharing
2. Extract the contents of the zip at _../Cell-Insights/_
3. Run main.py
4. Choose the image-sequences to run the piepline on. For sanity, choose the example .tif file from points 1 and 2 here above.
5. The piepline will ask for basic running parameters, the default is 'Yes' to all
   6. The trained models are _Incucyte_ and _INCell_, referring to the two imaging methods that were used in this study. The default usage is _INCell_
6. See the Tracker results at the chosen output directory. Default is _C:\CellInsights_

![](https://github.com/avivshalem/Cell-Insights/blob/main/cells.gif)

## Training your own Segmentation model
Please refer to the readme under the Segmentation_Pipeline directory

## The Feature Analysis
Please refer to the readme under the Feature_Pipeline directory

## The Classifier
Please refer to the readme under the Classifier_Pipeline directory

