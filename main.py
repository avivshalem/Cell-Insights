import os
import SegmentationPipeline as Segmentation
import TrackingPipeline as Tracking
import VisualizerExtractorPipeline as VisulizeExtractor
from tkinter import Tk, filedialog
import tifffile as tiff
import numpy as np
import time
import matplotlib.pyplot as plt


def main():
    print('Process Started.. \n')
    allTime = time.time()

    # Data selector
    root = Tk()
    root.withdraw()
    root.attributes('-topmost', True)
    filenames = filedialog.askopenfilenames(title='Open a TIF/TIFF file', filetypes=(('TIF files', '*.tif'), ('TIFF files', '*.tiff')))

    # Askers
    readerAsk = input("Run image separator? (YES/NO) ")
    segmentAsk = input("Run segmentation? (YES/NO) ")
    trackAsk = input("Run tracking? (YES/NO) ")
    extractAsk = input("Run extractor? (YES/NO) ")
    defaultAsk = input("Use default paths to save? \nNote: if NO, each experiment will need a specific path, and process will be stopped until path is chosen! (YES/NO) ")
    if defaultAsk == 'yes' or defaultAsk == 'Yes' or defaultAsk == 'YES':
        print("\nDefault path chosen, see outputs in C:\Tracking\Cell migration data\SegmentationAgain\n")

    sequenceOffset = {}
    zeroOffset = input("Use zero as offset for all experiments? (YES/NO) ")

    for filename in filenames:
        if not filename:
            print('Please choose a TIF/TIFF file!')
            exit()
        experimentName = filename[filename.rfind('/') + 1 : filename.rfind('.')]
        if zeroOffset:
            sequenceOffset[experimentName] = '0'
        else:
            sequenceOffset[experimentName] = input(f"Enter the sequence offset for experiment {experimentName}: ")

    for filename in filenames:
        experimentName = filename[filename.rfind('/') + 1 : filename.rfind('.')]
        print(f'\nExperiment {experimentName} Started..\n')
        experimentTime = time.time()

        # Folders selector
        if defaultAsk == 'yes' or defaultAsk == 'Yes' or defaultAsk == 'YES':
            saveFolder = ''
        else:
            saveFolder = filedialog.askdirectory(title='Choose a folder to save the Image Sequence')
        if not saveFolder:
            print('Image Sequence save folder was not chosen and default path is selected')
            saveFolder = r'C:\Tracking\Cell migration data\SegmentationAgain\Data'
        if experimentName not in saveFolder:
            saveFolderPNG = os.path.join(saveFolder, experimentName, 'PNG')
            saveFolder = os.path.join(saveFolder, experimentName, 'TIF')
            os.makedirs(saveFolder, exist_ok=True)
            os.makedirs(saveFolderPNG, exist_ok=True)

        if defaultAsk == 'yes' or defaultAsk == 'Yes' or defaultAsk == 'YES':
            resultFolder = ''
        else:
            resultFolder = filedialog.askdirectory(title='Choose a folder to save the Segmentation results')
        if not resultFolder:
            print('Segmentation results save folder was not chosen and default path is selected')
            resultFolder = r'C:\Tracking\Cell migration data\SegmentationAgain\Segmentation'
        if experimentName not in resultFolder:
            resultFolder = os.path.join(resultFolder, experimentName)
            os.makedirs(resultFolder, exist_ok=True)

        if defaultAsk == 'yes' or defaultAsk == 'Yes' or defaultAsk == 'YES':
            trackingFolder = ''
        else:
            trackingFolder = filedialog.askdirectory(title='Choose a folder to save the Tracking results')
        if not trackingFolder:
            print('Tracking results save folder was not chosen and default path is selected')
            trackingFolder = r'C:\Tracking\Cell migration data\SegmentationAgain\Tracking'
        if experimentName not in trackingFolder:
            trackingFolder = os.path.join(trackingFolder, experimentName)
            os.makedirs(trackingFolder, exist_ok=True)

        # Image reader and sequence writer
        if readerAsk == 'yes' or readerAsk == 'Yes' or readerAsk == 'YES':
            print('\nImage Separator Started..\n')
            seperateTime = time.time()
            imageSequence = tiff.imread(filename)
            for image in range(np.shape(imageSequence)[0]):
                imageIndex = image - int(sequenceOffset[experimentName])
                if imageIndex < 0:
                    imageIndex += np.shape(imageSequence)[0]
                tiff.imsave(os.path.join(saveFolder, f't{str(imageIndex).zfill(3)}.tif') , imageSequence[image])
                plt.imsave(os.path.join(saveFolderPNG, f't{str(imageIndex).zfill(3)}.png'), imageSequence[image], cmap='gray')
            print(f'Image Separator Done. Took {round((time.time() - seperateTime)/60, 2)} minutes\n')


        # Run segmentation
        if segmentAsk == 'yes' or segmentAsk == 'Yes' or segmentAsk == 'YES':
            print('\nSegmentation Started..\n')
            segmentTime = time.time()
            Segmentation.main(th_seed=7.2, th_cell=1.28, apply_clahe=True, savePath=saveFolder, resultPath=resultFolder, model_type='Incell', cuda=True, batch_size=1) #incell
            # Segmentation.main(th_seed=0.9*85, th_cell=0.08*70, apply_clahe=True, savePath=saveFolder, resultPath=resultFolder, model_type='Incocite', cuda=True, batch_size=1) #incocite
            print(f'\nSegmentation Done. Took {round((time.time() - segmentTime)/60, 2)} minutes\n')

        # Run tracking
        if trackAsk == 'yes' or trackAsk == 'Yes' or trackAsk == 'YES':
            print('\nTracking Started..\n')
            trackingTime = time.time()
            Tracking.main(img_path=saveFolder, segm_path=resultFolder, res_path=trackingFolder, delta_t=3, default_roi_size=5)
            print(f'\nTracking Done. Took {round((time.time() - trackingTime)/60, 2)} minutes\n')

        # Visualize tracks and Extract XY positions
        if extractAsk == 'yes' or extractAsk == 'Yes' or extractAsk == 'YES':
            VisulizeExtractor.main(trackingDir=trackingFolder, dataPath=saveFolderPNG)

        print(f'\nExperiment {experimentName} Done. Took {round((time.time() - experimentTime) / 60, 2)} minutes\n')

    # Finished
    print(f'\nEntire Process Completed. Took {round((time.time() - allTime)/60, 2)} minutes to process {len(filenames)} experiments\n')


if __name__ == "__main__":

    main()
