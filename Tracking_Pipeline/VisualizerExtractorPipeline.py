import os
import pickle
import numpy as np
import random
from moviepy.editor import ImageSequenceClip
import colorsys
from PIL import Image
import pandas as pd

def savePathImages(trackIndex, tracks, tracksAndTimes, traversed):
    traversed[0] += 1
    if traversed[0] > 1000000:
        print('traversed too long, breaking. ', trackIndex)
        return
    currentTrack = tracks[trackIndex]

    tracksAndTimes[trackIndex] = {**tracksAndTimes[trackIndex], **currentTrack.center_of_mass}


    if len(currentTrack.successors) == 0:
        return
    for nextTrackIndex in currentTrack.successors:
        try:
            if tracksAndTimes[nextTrackIndex]:
                continue
        except KeyError:
            tracksAndTimes[nextTrackIndex] = tracksAndTimes[trackIndex].copy()
            savePathImages(nextTrackIndex, tracks, tracksAndTimes, traversed)
    del tracksAndTimes[trackIndex]
    return


def main(trackingDir, dataPath):

    trackingDir = trackingDir
    dataPath = dataPath

    with open(os.path.join(trackingDir, f'tracks.pickle'), 'rb') as handle:
        tracks = pickle.load(handle)

    tracksAndTimes = {}
    traversed = {}
    for trackIndex in tracks:

        if tracks[trackIndex].pred_track_id != [0]:
            continue
        tracksAndTimes[trackIndex] = {}

        traversed[0] = -1
        savePathImages(trackIndex, tracks, tracksAndTimes, traversed)

    numColors = len(tracksAndTimes)
    HSV_tuples = [(x/numColors, 0.5, 0.5) for x in range(numColors)]
    RGB_tuples = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
    for rgb in range(len(RGB_tuples)):
        RGB_tuples[rgb] = tuple(255.0 * x for x in RGB_tuples[rgb]) + (0, )
    random.shuffle(RGB_tuples)

    temp = Image.open(os.path.join(dataPath, os.listdir(dataPath)[0]))
    img_dict = np.zeros((len(os.listdir(dataPath)),temp.size[1],temp.size[0],4), dtype=int)
    for i, item in enumerate(os.listdir(dataPath)):
        img_dict[i] = np.asarray(Image.open(os.path.join(dataPath, item)))


    for i in tracksAndTimes:
        if len(tracksAndTimes[i]) < 64:
            continue
        currC = RGB_tuples.pop()
        for t in tracksAndTimes[i]:
            img_dict[t, int(round(tracksAndTimes[i][t][0]) - 10):int(round(tracksAndTimes[i][t][0]) + 10),
            int(round(tracksAndTimes[i][t][1]) - 10):int(round(tracksAndTimes[i][t][1]) + 10)] = currC

    clip = ImageSequenceClip(list(img_dict), fps=5)
    clip.write_gif(os.path.join(trackingDir, 'gif.gif'), fps=5)

    df = pd.DataFrame.from_dict(tracksAndTimes)
    df = df.sort_index()
    df.to_excel(os.path.join(trackingDir, 'XY_Positions.xlsx'))

if __name__ == "__main__":

    main()