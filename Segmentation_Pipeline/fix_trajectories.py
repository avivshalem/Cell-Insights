import os

import pandas as pd
import numpy as np


def main(path):
    # path = r'C:\CellInsights\Tracking'
    pathList = os.listdir(path)
    dfNew = pd.DataFrame()
    logicalCounter = 0
    for experiment in pathList:
        if '.xlsx' in experiment:
            continue
        if '.csv' in experiment:
            continue
        try:
            df = pd.read_excel(os.path.join(path, experiment, r'XY_Positions.xlsx'), engine='openpyxl')
        except FileNotFoundError:
            continue    # if 'LatB' in experiment:
        #     HDorHC = 'HD_LatB'
        #     continue
        # elif 'DMSO' in experiment:
        #     HDorHC = 'HD_DMSO'
        #     continue
        # elif 'Media' in experiment:
        #     HDorHC = 'HD_Media'
        #     continue
        # elif 'GM' in experiment:
        #     HDorHC = 'HD'

            # if '2147' in experiment or '2147' in experiment or '4200' in experiment or '4196' in experiment or '4212' in experiment or '5031' in experiment or '6274' in experiment: # For Onset
            #     continue
        # elif experiment[0] == 'P':
        #     HDorHC = 'Progeria'
        #
        # elif 'fld' in experiment:
        #     HDorHC = 'LatB'
        #     continue
        # else:
        #     HDorHC = 'HC'
        #     continue
        # if 'GSD' in experiment:
        #      HDorHC = 'GSD'
        if 'HC - ' in experiment:
            HDorHC = 'HC'
        elif 'HD - ' in experiment:
            HDorHC = 'HD'
        elif 'LatBHigh - ' in experiment:
            # pass
            HDorHC = 'LatBHigh'
        elif 'LatBLow - ' in experiment:
            # pass
            HDorHC = 'LatBLow'
        elif 'HCMedia - ' in experiment:
            HDorHC = 'HCMedia'
        elif 'HDMedia - ' in experiment:
            HDorHC = 'HDMedia'
            # print("Note ", experiment)
        else:
            print("Note ", experiment)
        # HDorHC = 'HGPS'

        for i in range(1, df.shape[1]):
            currentCol = df[df.columns[i]].copy()
            if currentCol.count() < 75:
                continue
            currentCol = currentCol.dropna()

            logicalCounter += 1
            flag = False
            logicalCounterCells = 0
            dfNewInner = pd.DataFrame()
            good = True
            for item in currentCol:
                logicalCounterCells += 1
                item = [float(item.split(', ')[0].split('(')[1]), float(item.split(', ')[1].split(')')[0])]
                if flag:
                    distance = np.sqrt((item[0] - prevItem[0])**2 + (item[1] - prevItem[1])**2)
                    # print(distance)
                    if 'IC' in experiment: # Incell
                        dis = 150
                    else: # Incocite
                        dis = 150
                    if distance > dis:
                        print(i, ' has failed due to large movement')
                        logicalCounter -= 1
                        good = False
                        break
                if not flag:
                    flag = True
                prevItem = item.copy()
                a_row = pd.Series([HDorHC, experiment, item[0], item[1], logicalCounterCells, logicalCounter])
                row_df = pd.DataFrame([a_row])
                dfNewInner = pd.concat([dfNewInner, row_df], ignore_index=True)
            if good:
                dfNew = pd.concat([dfNew, dfNewInner], ignore_index=True)

    dfNew.to_csv(os.path.join(path, r'Output.csv'))


if __name__ == "__main__":
    main()