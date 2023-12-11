import os
import tifffile as tiff
import cv2
import numpy as np
import pandas as pd

filename = 'output_circularity'
path = r'C:\CellInsights\Segmentation'
path_list = os.listdir(path)
# dict_it = dict.fromkeys(['Experiment', dict.fromkeys(['Circularity', 'Circularity minimum enclosing circle', 'Circularity minimum enclosing ellipse',
#                          'Aspect Ratio', 'Eccentricity', 'Form Factor'])])
dict_it = {}
incell_res = 1/2000
incocite_res = 1/1500

res = incocite_res / incell_res

for folder in path_list:
    dict_it[folder] = {'Aspect Ratio': dict.fromkeys(['mean', 'median', 'std']),
                       'Eccentricity': dict.fromkeys(['mean', 'median', 'std']),
                       'Form Factor': dict.fromkeys(['mean', 'median', 'std'])}
    # dict_it[folder] = dict.fromkeys(['Circularity', 'Circularity minimum enclosing circle',
                                     # 'Circularity minimum enclosing ellipse', 'Aspect Ratio', 'Eccentricity',
                                     # 'Form Factor'])
    current_folder = os.path.join(path, folder)
    folder_list = os.listdir(current_folder)
    c1Mean, c2Mean, c3Mean = [], [], []
    c1Median, c2Median, c3Median = [], [], []
    c1STD, c2STD, c3STD= [], [], []
    for j, image in enumerate(folder_list):
        if j == 164:
            break
        current_image = os.path.join(current_folder, image)
        read_image = tiff.imread(current_image)
        read_image_normed = np.round((read_image/np.max(read_image)) * 255).astype(np.uint8)
        _, thresh = cv2.threshold(read_image_normed, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        contours, _ = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        c1, c2, c3 = [], [], []
        reduce3, reduce5 = 0, 0
        for i, contourOG in enumerate(contours):
            contour = (contourOG*res).astype(np.int32)
            area = cv2.contourArea(contour)
            perimeter = cv2.arcLength(contour, True)
            circularity = 4 * np.pi * (area / (perimeter * perimeter))
        #     c1.append(circularity)
        #     # print("Circularity:", circularity)

            (x, y), radius = cv2.minEnclosingCircle(contour)
            center = (int(x), int(y))
            radius = int(radius)
            A = cv2.contourArea(contour)
            circularity = A / (np.pi * radius * radius)
            # c2.append(circularity)
            # print("Circularity minimum enclosing circle:", circularity)

            if (len(contour) > 4):
                ellipse = cv2.fitEllipse(contour)
                A = cv2.contourArea(contour)
                a = ellipse[1][0] / 2
                b = ellipse[1][1] / 2
                circularity = A / (np.pi * a * b)
                # c3.append(circularity)
            else:
                reduce3 += 1
            # print("Circularity minimum enclosing ellipse:", circularity)

            rect = cv2.minAreaRect(contour)
            (x, y), (w, h), angle = rect

            # Calculate the aspect ratio
            aspect_ratio = min(w, h) / max(w, h)
            c1.append(aspect_ratio)
            # print("Aspect Ratio:", aspect_ratio)

            # Calculate the minimum bounding ellipse of the contour
            if (len(contour) > 4):
                ellipse = cv2.fitEllipse(contour)

                # Calculate the eccentricity
                a = max(ellipse[1][0], ellipse[1][1]) / 2
                b = min(ellipse[1][0], ellipse[1][1]) / 2
                eccentricity = np.sqrt(1 - (b * b) / (a * a))


                c2.append(eccentricity)
            else:
                reduce5 += 1
            # print("Eccentricity:", eccentricity)

            area = cv2.contourArea(contour)
            perimeter = cv2.arcLength(contour, True)
            # form_factor = (perimeter * perimeter) / area
            form_factor = (4*np.pi*area) / (perimeter*perimeter)
            c3.append(form_factor)
            # print("Form Factor:", form_factor)
        c1Mean.append(np.mean(c1))
        c2Mean.append(np.mean(c2))
        c3Mean.append(np.mean(c3))

        c1Median.append(np.median(c1))
        c2Median.append(np.median(c2))
        c3Median.append(np.median(c3))

        c1STD.append(np.std(c1))
        c2STD.append(np.std(c2))
        c3STD.append(np.std(c3))


    ck1 = [np.mean(c1Mean), np.mean(c2Mean), np.mean(c3Mean)]
    ck2 = [np.median(c1Median), np.median(c2Median), np.median(c3Median)]
    ck3 = [np.std(c1STD), np.std(c2STD), np.std(c3STD)]

    for k, key in enumerate(dict_it[folder].keys()):
        dict_it[folder][key]['mean'] = ck1[k]
        dict_it[folder][key]['median'] = ck2[k]
        dict_it[folder][key]['std'] = ck3[k]
df = pd.DataFrame.from_dict(dict_it, orient='index')
df.to_excel(path + f'\{filename}.xlsx')

