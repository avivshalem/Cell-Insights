import os

import pandas as pd
import matplotlib.pyplot as plt
import json
import numpy as np
import scipy.stats as stats
import natsort
from scipy.stats import shapiro, kstest, anderson, norm
from statsmodels.graphics.gofplots import qqplot
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from collections import defaultdict
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout
from sklearn.preprocessing import LabelEncoder
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from itertools import chain, combinations
from keras.callbacks import ModelCheckpoint
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.utils import resample
from scipy import interp
import math

def is_close(a, b, tolerance=1e-3):
    return abs(a - b) > tolerance

def flatten_dict(d):
    result = []
    for key in d:
        result.extend(d[key])
    return result

def check_normal_dist(dic):
    # Perform the Shapiro-Wilk test
    data = np.asarray(flatten_dict(dic))

    # Visualize the data
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].hist(data, bins=30)
    axes[0].set_title("Histogram")

    axes[1].plot(np.sort(data), norm.cdf(np.sort(data)))
    axes[1].set_title("Normal Probability Plot")

    qqplot(data, line='s', ax=axes[2])
    axes[2].set_title("Q-Q Plot")

    stat, p = shapiro(data)
    alpha = 0.05
    res = 0
    # print(f"Shapiro-Wilk test: stat={stat:.3f}, p-value={p:.3f}")
    if p > alpha:
        res += 1
        # print("Data is normally distributed")
    else:
        pass
        # print("Data is not normally distributed")

    # # Perform the Kolmogorov-Smirnov test
    # stat, p = kstest(data, 'norm')
    # # print(f"Kolmogorov-Smirnov test: stat={stat:.3f}, p-value={p:.3f}")
    # if p > alpha:
    #     res += 1
    #     # print("Data is normally distributed")
    # else:
    #     pass
        # print("Data is not normally distributed")

    # result = anderson(data)
    #
    # # print("Anderson-Darling test results:")
    # # print(f"Statistic: {result.statistic:.3f}")
    # # print(f"Critical values: {result.critical_values}")
    # # print(f"Significance level: {result.significance_level:.3f}")
    #
    # if result.statistic > result.critical_values[2]:
    #     pass
    #     # print("Data is not normally distributed")
    # else:
    #     res += 1
    #     # print("Data is normally distributed")

    if res == 1:
        print('Data is normally distributed with a decision of 1 test')
        return True
    if res == 2:
        print('Data is normally distributed with a decision of 2 tests')
        return True
    if res == 3:
        print('Data is normally distributed with a decision of all 3 tests')
        return True
    print('Data is not normally distributed')
    return False

def remove_outliers(data):
    data = np.asarray(data)
    # model = IsolationForest(contamination=0.33)
    # model.fit(data.reshape(-1, 1))
    # outliers = model.predict(data.reshape(-1, 1))
    # data_clean = data[outliers == 1]
    # return data_clean
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    data = np.array(data)
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    mask = (data >= lower_bound) & (data <= upper_bound)
    return data[mask]

def extract_pValue(dic):
    normal = check_normal_dist(dic)
    for key1 in dic.keys():
        for key2 in dic.keys():
            if key1 != key2:
                # Define your null and alternative hypotheses
                null_hypothesis = "There is no difference between the means of {} and {}".format(key1, key2)
                alternative_hypothesis = "There is a significant difference between the means of {} and {}".format(key1,
                                                                                                                   key2)
                # Conduct a t-test
                t_statistic, p_value = stats.ttest_ind(dic[key1], dic[key2])
                u_statistic, p_valueMW = stats.mannwhitneyu(dic[key1], dic[key2])

                # Print the results
                if p_value < 0.05 and normal:
                    print("{} vs {}: t-statistic = {:.2f}, p-value = {:.4f}".format(key1, key2, t_statistic, p_value))
                if p_valueMW < 0.05 and not normal:
                    print("{} vs {}: U-statistic = {:.2f}, p-value = {:.4f}".format(key1, key2, u_statistic, p_valueMW))

# Load the Excel file into a Pandas DataFrame
file = 'PlateResults - LaminB+a actinin- S exp4 05.07.22__2022-07-12T10_09_44-Measurement 1'
df1 = pd.read_excel(rf'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\LaminB Data\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- S exp4 05.07.22 (9-10)__2022-07-12T11_44_17-Measurement 1'
df2 = pd.read_excel(rf'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\LaminB Data\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- R exp4 05.07.22__2022-07-17T09_42_44-Measurement 1'
df3 = pd.read_excel(rf'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\LaminB Data\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- exp8 12.02.23__2023-02-20T14_52_23-Measurement 1'
df4 = pd.read_excel(rf'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\LaminB Data\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- exp7 06.02.23__2023-02-19T14_31_42-Measurement 1'
df5 = pd.read_excel(rf'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\LaminB Data\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- exp6  05.02.23__2023-02-19T12_14_21-Measurement 1'
df6 = pd.read_excel(rf'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\LaminB Data\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- exp5 plate1 08.01.23__2023-01-12T15_50_14-Measurement 1'
df7 = pd.read_excel(rf'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\LaminB Data\{file}.xlsx')
file = 'HGPS PlateResults'
df8 = pd.read_excel(rf'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\LaminB Data\{file}.xlsx')
dfs = [df1, df3, df4, df5, df6, df7, df8]
dfsNew = []
dfCached = pd.concat(dfs, axis=0)
# unqiueCAPS = []
counter = 0
for i, d in enumerate(dfs):
    print('***********************************************')
    print(i)
    print('***********************************************')
    print(np.unique(d['CAP']))
    # if (i == 3 or i == 5 or i == 6): # filter CAP 0
    #     d = d[is_close(d['CAP'], 0)]
    # if (i == 4): # filter CAP 60
    #     d = d[is_close(d['CAP'], 60.24653313)]
    # if (i == 5): # filter CAP 69
    #     d = d[is_close(d['CAP'], 69.33744222)]
    # if ( i == 5 or i == 6): # filter CAP 74
    #     d = d[is_close(d['CAP'], 74.57627119)]
    # if ( i == 4): # filter CAP 76
    #     d = d[is_close(d['CAP'], 76.27118644)]
    # if (i == 0):  # filter CAP 126
    #     d = d[is_close(d['CAP'], 126.34822804)]
    # if (i == 2):  # filter CAP 129
    #     d = d[is_close(d['CAP'], 129.42989214)]
    # if (i == 2):  # filter CAP 140
    #     d = d[is_close(d['CAP'], 140.5238829)]
    # if (i == 6):  # filter CAP 140
    #     d = d[is_close(d['CAP'], 84.12942989)]
    #     d = d[is_close(d['CAP'], 96.14791988)]
    # if (i == 0):  # filter CAP 140
    #     d = d[is_close(d['CAP'], 114.02157165)]

    # d['CAP'] = (d['CAP']*1000)+counter
    counter += 10

    dfsNew.append(d)
    # if i == 0:
    #     dfCached['CAP'] = (dfCached['CAP']*1000)+counter
    # else:
    #     dfCached['CAP'] = (dfCached['CAP'])+counter
    # currentCAP = np.unique(dfCached['CAP'].tolist()).tolist()
    # unqiueCAPS.extend(currentCAP)
# unqiueCAPS = np.unique(unqiueCAPS).tolist()
# print(sorted(unqiueCAPS))

df = pd.concat(dfsNew, axis=0)
features = df.columns[4:-6].to_list()
# counter = 0
feature_names = []
for fi, feature in enumerate(features):
    # if "chanel_4" not in feature:
    if "chanel_4" not in feature and "chanel_1 Ratio Width to Length" not in feature:

        continue
    feature_names.append(feature)
    # if "GABOR" in feature or "SER" in feature or "Gabor" in feature:
    #     continue
    try:
        norm_to = df1[feature].mean()
    except KeyError:
        continue
    grouped = df1.groupby('CAP')[feature].apply(lambda x: x.values.tolist())
    new_min, new_max = min(grouped[0]), max(grouped[0])

    # grouped2 = df2.groupby('CAP')[feature].apply(lambda x: x.values.tolist())
    # grouped2Saved = grouped2.copy()
    # data_min, data_max = min(grouped2[0]), max(grouped2[0])
    # scaling_factor = (new_max - new_min) / (data_max - data_min)
    # shift = new_min - data_min * scaling_factor
    # for key, values in grouped2.items():
    #     grouped2[key] = np.array(values) * scaling_factor + shift

    df3Old = df3.copy()
    grouped3 = df3.groupby('CAP')[feature].apply(lambda x: x.values.tolist())
    grouped3Saved = grouped3.copy()
    data_min, data_max = min(grouped3[0]), max(grouped3[0])
    scaling_factor = (new_max - new_min) / (data_max - data_min)
    shift = new_min - data_min * scaling_factor
    for key, values in grouped3.items():
        grouped3[key] = np.array(values) * scaling_factor + shift
    temp_df = pd.DataFrame()
    for key in grouped3.keys():
        indices = df3[df3['CAP'] == key].index
        temp_df = temp_df.append(pd.DataFrame({feature: grouped3[key], 'index': indices}), ignore_index=True)
    feature_column_location = df3.columns.get_loc(feature)
    df3 = df3.drop(columns=[feature])
    df3 = df3.reset_index(drop=True)
    temp_df = temp_df.sort_values('index').reset_index(drop=True)
    df3.insert(feature_column_location, feature, temp_df[feature])

    grouped4 = df4.groupby('CAP')[feature].apply(lambda x: x.values.tolist())
    grouped4Saved = grouped4.copy()
    data_min, data_max = min(grouped4[0]), max(grouped4[0])
    scaling_factor = (new_max - new_min) / (data_max - data_min)
    shift = new_min - data_min * scaling_factor
    for key, values in grouped4.items():
        grouped4[key] = np.array(values) * scaling_factor + shift
    temp_df = pd.DataFrame()
    for key in grouped4.keys():
        indices = df4[df4['CAP'] == key].index
        temp_df = temp_df.append(pd.DataFrame({feature: grouped4[key], 'index': indices}), ignore_index=True)
    feature_column_location = df4.columns.get_loc(feature)
    df4 = df4.drop(columns=[feature])
    df4 = df4.reset_index(drop=True)
    temp_df = temp_df.sort_values('index').reset_index(drop=True)
    df4.insert(feature_column_location, feature, temp_df[feature])

    grouped5 = df5.groupby('CAP')[feature].apply(lambda x: x.values.tolist())
    grouped5Saved = grouped5.copy()
    data_min, data_max = min(grouped5[0]), max(grouped5[0])
    scaling_factor = (new_max - new_min) / (data_max - data_min)
    shift = new_min - data_min * scaling_factor
    for key, values in grouped5.items():
        grouped5[key] = np.array(values) * scaling_factor + shift
    temp_df = pd.DataFrame()
    for key in grouped5.keys():
        indices = df5[df5['CAP'] == key].index
        temp_df = temp_df.append(pd.DataFrame({feature: grouped5[key], 'index': indices}), ignore_index=True)
    feature_column_location = df5.columns.get_loc(feature)
    df5 = df5.drop(columns=[feature])
    df5 = df5.reset_index(drop=True)
    temp_df = temp_df.sort_values('index').reset_index(drop=True)
    df5.insert(feature_column_location, feature, temp_df[feature])

    grouped6 = df6.groupby('CAP')[feature].apply(lambda x: x.values.tolist())
    grouped6Saved = grouped6.copy()
    data_min, data_max = min(grouped6[0]), max(grouped6[0])
    scaling_factor = (new_max - new_min) / (data_max - data_min)
    shift = new_min - data_min * scaling_factor
    for key, values in grouped6.items():
        grouped6[key] = np.array(values) * scaling_factor + shift
    temp_df = pd.DataFrame()
    for key in grouped6.keys():
        indices = df6[df6['CAP'] == key].index
        temp_df = temp_df.append(pd.DataFrame({feature: grouped6[key], 'index': indices}), ignore_index=True)
    feature_column_location = df6.columns.get_loc(feature)
    df6 = df6.drop(columns=[feature])
    df6 = df6.reset_index(drop=True)
    temp_df = temp_df.sort_values('index').reset_index(drop=True)
    df6.insert(feature_column_location, feature, temp_df[feature])

    grouped7 = df7.groupby('CAP')[feature].apply(lambda x: x.values.tolist())
    grouped7Saved = grouped7.copy()
    data_min, data_max = min(grouped7[0]), max(grouped7[0])
    scaling_factor = (new_max - new_min) / (data_max - data_min)
    shift = new_min - data_min * scaling_factor
    for key, values in grouped7.items():
        grouped7[key] = np.array(values) * scaling_factor + shift
    temp_df = pd.DataFrame()
    for key in grouped7.keys():
        indices = df7[df7['CAP'] == key].index
        temp_df = temp_df.append(pd.DataFrame({feature: grouped7[key], 'index': indices}), ignore_index=True)
    feature_column_location = df7.columns.get_loc(feature)
    df7 = df7.drop(columns=[feature])
    df7 = df7.reset_index(drop=True)
    temp_df = temp_df.sort_values('index').reset_index(drop=True)
    df7.insert(feature_column_location, feature, temp_df[feature])

    grouped8 = df8.groupby('CAP')[feature].apply(lambda x: x.values.tolist())
    grouped8Saved = grouped8.copy()
    data_min, data_max = min(grouped8[0]), max(grouped8[0])
    scaling_factor = (new_max - new_min) / (data_max - data_min)
    shift = new_min - data_min * scaling_factor
    for key, values in grouped8.items():
        grouped8[key] = np.array(values) * scaling_factor + shift
    temp_df = pd.DataFrame()
    for key in grouped8.keys():
        indices = df8[df8['CAP'] == key].index
        temp_df = temp_df.append(pd.DataFrame({feature: grouped8[key], 'index': indices}), ignore_index=True)
    feature_column_location = df8.columns.get_loc(feature)
    df8 = df8.drop(columns=[feature])
    df8 = df8.reset_index(drop=True)
    temp_df = temp_df.sort_values('index').reset_index(drop=True)
    df8.insert(feature_column_location, feature, temp_df[feature])

dfs = [df1, df3, df4, df5, df6, df7, df8]
dfsNew = []
dfCached = pd.concat(dfs, axis=0)
# unqiueCAPS = []
counter = 0
for i, d in enumerate(dfs):
    print('***********************************************')
    print(i)
    print('***********************************************')
    print(np.unique(d['CAP']))
    # if (i == 3 or i == 5 or i == 6): # filter CAP 0
    #     d = d[is_close(d['CAP'], 0)]
    # if (i == 4): # filter CAP 60
    #     d = d[is_close(d['CAP'], 60.24653313)]
    # if (i == 5): # filter CAP 69
    #     d = d[is_close(d['CAP'], 69.33744222)]
    # if ( i == 5 or i == 6): # filter CAP 74
    #     d = d[is_close(d['CAP'], 74.57627119)]
    # if ( i == 4): # filter CAP 76
    #     d = d[is_close(d['CAP'], 76.27118644)]
    # if (i == 0):  # filter CAP 126
    #     d = d[is_close(d['CAP'], 126.34822804)]
    # if (i == 2):  # filter CAP 129
    #     d = d[is_close(d['CAP'], 129.42989214)]
    # if (i == 2):  # filter CAP 140
    #     d = d[is_close(d['CAP'], 140.5238829)]
    # if (i == 6):  # filter CAP 140
    #     d = d[is_close(d['CAP'], 84.12942989)]
    #     d = d[is_close(d['CAP'], 96.14791988)]
    # if (i == 0):  # filter CAP 140
    #     d = d[is_close(d['CAP'], 114.02157165)]

    # d['CAP'] = (d['CAP']*1000)+counter
    counter += 10

    dfsNew.append(d)
    # if i == 0:
    #     dfCached['CAP'] = (dfCached['CAP']*1000)+counter
    # else:
    #     dfCached['CAP'] = (dfCached['CAP'])+counter
    # currentCAP = np.unique(dfCached['CAP'].tolist()).tolist()
    # unqiueCAPS.extend(currentCAP)
# unqiueCAPS = np.unique(unqiueCAPS).tolist()
# print(sorted(unqiueCAPS))

df = pd.concat(dfsNew, axis=0)
features = df.columns[4:-6].to_list()
# counter = 0

# Group the data by the values in column A and aggregate the other columns by taking their mean
medianprops = dict(linestyle='-', linewidth=1.5, color='black')

features = df.columns[4:-6].to_list()
# counter = 0
dic_lists = []
value_lists = []
for fi, feature in enumerate(features):
    if "chanel_4" not in feature and "chanel_1 Ratio Width to Length" not in feature:
    # if "chanel_4" not in feature:
        continue
    # if "GABOR" in feature or "SER" in feature or "Gabor" in feature:
    #     continue
    # if "chanel_1 Ratio Width to Length" in feature:
    # if "spots_chanel_3_final - Spot Area [px²] - Mean per Well" not in feature:
    #     continue
    # try:
    #     norm_to = df7[feature].mean()
    # except KeyError:
    #     continue
    # if 'ir_chanel_4_at - ir_chanel_4_at Area' not in feature:
    #     continue
    # counter += 1
    print(f'***********************{feature}***********************')
    grouped = df.groupby('CAP')[feature].apply(lambda x: x.values.tolist())
    # grouped = grouped.apply(lambda x: [n / norm_to for n in x])
    groupedRaw = grouped.copy()
    groupedOutliers = grouped.copy()
    groupedLog = grouped.copy()
    HCItems = []
    PItems = []
    MItems = []
    SItems = []
    HGPSItems = []
    HCItemsRaw = []
    PItemsRaw = []
    MItemsRaw = []
    SItemsRaw = []
    HGPSItemsRaw = []
    HCItemsOutliers = []
    PItemsOutliers = []
    MItemsOutliers = []
    SItemsOutliers = []
    HGPSItemsOutliers = []
    HCItemsLog = []
    PItemsLog = []
    MItemsLog = []
    SItemsLog = []
    HGPSItemsLog = []
    for cap in grouped.keys():
        # beforeLength = len(grouped[cap])
        groupedRaw[cap] = grouped[cap]
        groupedOutliers[cap] = remove_outliers(grouped[cap])
        grouped[cap] = remove_outliers(grouped[cap])
        groupedLog[cap] = np.log(grouped[cap])
        grouped[cap] = grouped[cap]
        # afterLength = len(grouped[cap])
        # if beforeLength != afterLength:
        #     print('Outlier removed for cap score ', cap)
        if cap == 0:
            HCItemsRaw += list(groupedRaw[cap])
            HCItemsOutliers += list(groupedOutliers[cap])
            HCItemsLog += list(groupedLog[cap])
            HCItems += list(grouped[cap])
        elif 0 < cap < 90:
            PItemsRaw += list(groupedRaw[cap])
            PItemsOutliers += list(groupedOutliers[cap])
            PItemsLog += list(groupedLog[cap])
            PItems += list(grouped[cap])
        elif 90 <= cap < 114:
            MItemsRaw += list(groupedRaw[cap])
            MItemsOutliers += list(groupedOutliers[cap])
            MItemsLog += list(groupedLog[cap])
            MItems += list(grouped[cap])
        elif cap == 999:
            HGPSItemsRaw += list(groupedRaw[cap])
            HGPSItemsOutliers += list(groupedOutliers[cap])
            HGPSItemsLog += list(groupedLog[cap])
            HGPSItems += list(grouped[cap])
        else:
            SItemsRaw += list(groupedRaw[cap])
            SItemsOutliers += list(groupedOutliers[cap])
            SItemsLog += list(groupedLog[cap])
            SItems += list(grouped[cap])
    # HCItems = remove_outliers(HCItems)
    # PItems = remove_outliers(PItems)
    # MItems = remove_outliers(MItems)
    # SItems = remove_outliers(SItems)

    mean = np.mean(HCItems)
    std = np.std(HCItems)
    dicRaw = {'HC': HCItemsRaw, "Premanifest": PItemsRaw, "Mild": MItemsRaw, "Severe": SItemsRaw, "HGPS": HGPSItemsRaw}
    dicOutliers = {'HC': HCItemsOutliers, "Premanifest": PItemsOutliers, "Mild": MItemsOutliers, "Severe": SItemsOutliers, "HGPS": HGPSItemsOutliers}
    dicLog = {'HC': HCItemsLog, "Premanifest": PItemsLog, "Mild": MItemsLog, "Severe": SItemsLog, "HGPS": HGPSItemsLog}
    # dic = {'HC': ((HCItems-mean)/std).tolist(), "Premanifest": ((PItems-mean)/std).tolist(), "Mild": ((MItems-mean)/std).tolist(), "Severe": ((SItems-mean)/std).tolist()}
    dic = {'HC': (HCItems / mean).tolist(), "Premanifest": (PItems / mean).tolist(), "Mild": (MItems / mean).tolist(), "Severe": (SItems / mean).tolist(), "HGPS": (HGPSItems / mean).tolist()}
    l = []
    l.extend((PItems / mean).tolist())
    l.extend((MItems / mean).tolist())
    l.extend((SItems / mean).tolist())
    # dic = {'HC': (HCItems / mean).tolist(), "HD": l, "HGPS": (HGPSItems / mean).tolist()}
    # l.extend(PItemsRaw)
    # l.extend(MItemsRaw)
    # l.extend(SItemsRaw)
    # dicRaw = {'HC': HCItemsRaw,"HD": l, "HGPS": HGPSItemsRaw}
    # l.extend(PItemsLog)
    # l.extend(MItemsLog)
    # l.extend(SItemsLog)
    # dicLog = {'HC': HCItemsLog,"HD": l, "HGPS": HGPSItemsLog}

    dic_lists.append(dic)
    grouped_value = (grouped.values) / mean
    value_lists.append(grouped_value)

def flatten_dict(d):
    result = []
    for key in d:
        result.extend(d[key])
    return result

def check_normal_dist(dic):
    # Perform the Shapiro-Wilk test
    data = np.asarray(flatten_dict(dic))

    # Visualize the data
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].hist(data, bins=30)
    axes[0].set_title("Histogram")

    axes[1].plot(np.sort(data), norm.cdf(np.sort(data)))
    axes[1].set_title("Normal Probability Plot")

    qqplot(data, line='s', ax=axes[2])
    axes[2].set_title("Q-Q Plot")

    # stat, p = shapiro(data)
    alpha = 0.05
    res = 0
    # # print(f"Shapiro-Wilk test: stat={stat:.3f}, p-value={p:.3f}")
    # if p > alpha:
    #     res += 1
    #     # print("Data is normally distributed")
    # else:
    #     pass
    #     # print("Data is not normally distributed")

    # Perform the Kolmogorov-Smirnov test
    # stat, p = kstest(data, 'norm')
    # # print(f"Kolmogorov-Smirnov test: stat={stat:.3f}, p-value={p:.3f}")
    # if p > alpha:
    #     res += 1
    #     # print("Data is normally distributed")
    # else:
    #     pass
    #     # print("Data is not normally distributed")
    #
    result = anderson(data)
    #
    # print("Anderson-Darling test results:")
    # print(f"Statistic: {result.statistic:.3f}")
    # print(f"Critical values: {result.critical_values}")
    # print(f"Significance level: {result.significance_level:.3f}")

    if result.statistic > result.critical_values[2]:
        pass
        # print("Data is not normally distributed")
    else:
        res += 1
        # print("Data is normally distributed")
    if res == 1:
        print('Data is normally distributed with a decision of 1 test')
        return True
    if res == 2:
        print('Data is normally distributed with a decision of 2 tests')
        return True
    if res == 3:
        print('Data is normally distributed with a decision of all 3 tests')
        return True
    print('Data is not normally distributed')
    return False

def alignment_values_to_all(dic, key_to_align):
    # Calculate the minimum and maximum values in the key to align to
    min_to_align_to = min(min(dic['Premanifest'], dic['Mild'], dic['Severe'], dic['HGPS']))
    max_to_align_to = max(max(dic['Premanifest'], dic['Mild'], dic['Severe'], dic['HGPS']))

    # Calculate the minimum and maximum values in the key to align
    min_to_align = min(dic[key_to_align])
    max_to_align = max(dic[key_to_align])

    # Calculate the scaling factor to align the values in the key to align to the range of values in the key to align
    scaling_factor = (max_to_align - min_to_align) / (max_to_align_to - min_to_align_to)

    # Loop through each value in the key to align
    for i, value in enumerate(dic[key_to_align]):
        # Scale the value based on the scaling factor
        dic[key_to_align][i] = (value - min_to_align) / scaling_factor + min_to_align_to

    return dic

def alignment_values(dic, key_to_align, key_to_align_to, lower=25, upper=75):
    # Calculate the minimum and maximum values in the key to align to
    min_to_align_to = min(dic[key_to_align_to])
    max_to_align_to = max(dic[key_to_align_to])

    # Calculate the minimum and maximum values in the key to align
    min_to_align = min(dic[key_to_align])
    max_to_align = max(dic[key_to_align])

    # Calculate the scaling factor to align the values in the key to align to the range of values in the key to align
    scaling_factor = (max_to_align - min_to_align) / (max_to_align_to - min_to_align_to)

    # Loop through each value in the key to align
    for i, value in enumerate(dic[key_to_align]):
        # Scale the value based on the scaling factor
        dic[key_to_align][i] = (value - min_to_align) / scaling_factor + min_to_align_to

    return dic
def extract_pValue(dic):
    normal = check_normal_dist(dic)
    for key1 in dic.keys():
        for key2 in dic.keys():
            if key1 != key2:
                # Define your null and alternative hypotheses
                null_hypothesis = "There is no difference between the means of {} and {}".format(key1, key2)
                alternative_hypothesis = "There is a significant difference between the means of {} and {}".format(key1,
                                                                                                                   key2)

                # Conduct a t-test
                t_statistic, p_value = stats.ttest_ind(dic[key1], dic[key2])
                u_statistic, p_valueMW = stats.mannwhitneyu(dic[key1], dic[key2])

                # Print the results

                if p_value < 0.05 and normal:
                    print("{} vs {}: t-statistic = {:.2f}, p-value = {:.4f}".format(key1, key2, t_statistic, p_value))
                if p_valueMW < 0.05 and not normal:
                    print("{} vs {}: U-statistic = {:.2f}, p-value = {:.4f}".format(key1, key2, u_statistic, p_valueMW))

CounterXL = 0
df1 = pd.read_excel(r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\output.xlsx')
df2 = pd.read_excel(r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\output707_708_24H_normed.xlsx')
df3 = pd.read_excel(r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\outputHCR.xlsx')
df = pd.concat([df1, df2, df3])
df_new = pd.DataFrame()
cols = df.iloc[:, 0]
new_cols = []
for itemIndex, item in enumerate(cols):
    if '180220cellmigprogeria' in item:
        if '- 02' in item:
            new_cols.append('999 - HGPS - pro271/13')
        if '- 03' in item:
            new_cols.append('999 - HGPS - pro367/13')
        if '- 04' in item:
            new_cols.append('999 - HGPS - pro169/15')
        if '- 07' in item:
            new_cols.append('999 - HGPS - pro127/16')
        if '- 09' in item:
            new_cols.append('999 - HGPS - pro122/18')
        if '- 05' in item:
            continue
            # new_cols.append('EV/16')
        if '- 06' in item:
            continue
            # new_cols.append('AV/17')
        if '- 08' in item:
            continue
            # new_cols.append('GSD1a53273/9')
        if '- 10' in item:
            continue
            # new_cols.append('Unknown')
    elif 'HCMIG21012020' in item:
        if '- 02' in item: # 0
            # continue
            new_cols.append('0 - HC - NA1016-7')
        if '- 03' in item: # 1
            continue
            new_cols.append('0 - HC - NA0878-9')
        if '- 04' in item: # 2
            continue
            new_cols.append('0 - HC - NA0025-8')
        if '- 07' in item: # 3
            continue
            new_cols.append('0 - HC - NA0233-8')
        if '- 09' in item: # 4
            continue
            new_cols.append('0 - HC - NA0561-8')
        if '- 05' in item: # 5
            continue
            new_cols.append('0 - HC - NA1170-5')
        if '- 06' in item: # 6
            continue
            new_cols.append('0 - HC - NA0633-6')
        if '- 08' in item: # 7
            continue
            new_cols.append('0 - HC - NA0579-10')
        if '- 10' in item: # 8
            continue
            new_cols.append('0 - HC - NA0769-7')
        if '- 11' in item: # 9
            continue
            new_cols.append('0 - HC - NA0143-10')

    elif 'HD migration' in item:
        if '- 02' in item:
            new_cols.append('111 - Mild - GM04887/8')
        if '- 03' in item:
            new_cols.append('100 - Mild - GM04212/8')
        if '- 04' in item:
            new_cols.append('132 - Severe - GM04476/10')
        if '- 07' in item:
            new_cols.append('96 - Mild - GM04819/8')
        if '- 09' in item:
            new_cols.append('110 - Mild - GM04196/8')
        if '- 05' in item:
            new_cols.append('129 - Severe - GM06274/8')
        if '- 06' in item:
            new_cols.append('119 - Severe - GM02147/12')
        if '- 08' in item:
            new_cols.append('94 - Mild - GM04799/7')
    elif 'VID707' in item:
        if '3_' in item: # 10
            # continue
            new_cols.append('0 - HC - NA1170')
        elif '4_' in item:
            new_cols.append('76 - Premanifest - GM04693')
        elif '5_' in item:
            if 'E' in item or 'F' in item or 'G' in item: # 11
                continue
                new_cols.append('0 - HC - AG16146')
            else: # 12
                # continue
                new_cols.append('0 - HC - GM01653')
        elif '6_' in item:
            new_cols.append('76 - Premanifest - GM04847')
        elif '2_' in item:
            new_cols.append('60 - Premanifest - GM04837')
        elif '8_' in item:
            new_cols.append('129 - Severe -GM00305')
        elif '9_' in item:
            if 'E' in item or 'F' in item or 'G' in item: # 13
                # continue
                new_cols.append('0 - HC - NA0971')
            else: # 14
                # continue
                new_cols.append('0 - HC - NA0495')
        elif '10_' in item:
            new_cols.append('141 - Severe - GM02165')
        else:  # 7
            new_cols.append('69 - Premanifest - GM04689')
    elif 'VID708' in item:
        if '3_' in item:
            if 'E' in item or 'F' in item or 'G' in item: # 15
                # continue
                new_cols.append('0 - HC - GM01653')
            else: # 16
                # continue
                new_cols.append('0 - HC - GM01650')
        elif '4_' in item:
            new_cols.append('115 - Severe - GM04691')
        elif '5_' in item: # 17
            # continue
            new_cols.append('0 - HC - NA1170')
        elif '6_' in item:
            new_cols.append('104 - Mild - GM04849')
        elif '7_' in item:
            new_cols.append('126 - Severe - GM04807')
        elif '2_' in item:
            new_cols.append('114 - Severe - GM04687')
        elif '9_' in item:
            new_cols.append('106 - Mild - GM04767')
        elif '10_' in item:
            new_cols.append('126 - Severe - GM04287')
        else:  # 18
            if 'E' in item or 'F' in item or 'G' in item: # 18
                # continue
                new_cols.append('0 - HC - NA0143')
            else: # 19
                # continue
                new_cols.append('0 - HC - NA0848')
    elif 'HC - NA' in item:
        # continue
        new_cols.append("0 - " + item)
        # new_cols.append(str(25+CounterXL) + " - " + item)
        CounterXL += 1
    else:
        if '- 08' in item:
            new_cols.append('131 - Severe - GM04200')
        if '- 02' in item:
            if 'E - ' in item or 'F - ' in item or 'G - ' in item: # 20
                # continue
                new_cols.append('0 - HC - GM04287F2.1')
            else: # 21
                # continue
                new_cols.append('0 - HC - GM04709B2.1')
        if '- 05' in item:
            new_cols.append('111 - Mild - GM04715')
        if '- 03' in item:
            new_cols.append('75 - Premanifest - GM04717')
        if '- 11' in item:
            new_cols.append('84 - Premanifest - GM04719')
        if '- 06' in item:
            new_cols.append('91 - Mild - GM04721')
        if '- 09' in item:
            new_cols.append('139 - Severe - GM05031')
        if '- 07' in item: # 22
            # continue
            new_cols.append('0 - HC - NA0157')
        if '- 10' in item: # 23
            # continue
            new_cols.append('0 - HC - NA0234')
        if '- 04' in item: # 24
            # continue
            new_cols.append('0 - HC - NA0477')
    df_new = df_new.append(df.iloc[itemIndex, :])

type_cols = []
for col in range(len(new_cols)):
    if 'HC' in new_cols[col]:
        type_cols.append('HC')
    elif 'Premanifest' in new_cols[col]:
        type_cols.append('Premanifest')
    elif 'Mild' in new_cols[col]:
        type_cols.append('Mild')
    elif 'Severe' in new_cols[col]:
        type_cols.append('Severe')
    elif 'HGPS' in new_cols[col]:
        type_cols.append('HGPS')
    new_cols[col] = new_cols[col][:new_cols[col].find(" ")]


df_new = df_new.iloc[:, 1:]
df1 = df_new
df1.insert(0,'Unnamed: 0',type_cols)
df1.set_index(df1.columns[0], inplace=True)
df1 = df1.applymap(lambda x: x.replace("'", '"') if isinstance(x, str) else x)

# Convert the strings to dictionaries
df1 = df1.applymap(lambda x: json.loads(x) if isinstance(x, str) else x if pd.notnull(x) else np.nan)
mean = {}
std = {}
all = {}
allNone = {}
allCached = {}
allLog = {}
for item in df1.columns:
    dicMean = {}
    dicMeanNone = {}
    dicMeanLog = {}
    dicMedian = {}
    dicStd = {}
    # fig, ax = plt.subplots()
    a = df1[item]
    for index, row in a.items():
        for stat in ['mean', 'median', 'std']:
            if stat == 'mean':
                if index not in dicMean:
                    dicMean[index] = [row[stat]]
                else:
                    dicMean[index].append(row[stat])

    dicMeanCached = dicMean.copy()
    for key in dicMean:
        data = np.asarray(dicMean[key])

        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        data = np.array(data)
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        mask = (data >= lower_bound) & (data <= upper_bound)
        dicMean[key] = data[mask].tolist()

    for key in dicMean:
        dicMeanNone[key] = dicMean[key]
        dicMeanLog[key] = np.log(dicMean[key])
        # dicMean[key] = np.log(dicMean[key])
    # extract_pValue(dicStd)
    mean[item] = np.mean(dicMean['HC'])
    std[item] = np.std(dicMean['HC'])
    for key in dicMean:
        # dicMean[key] = (dicMean[key] - mean[item]) / std[item]
        dicMean[key] = dicMean[key] / mean[item]
    print(f'***************************{item}************************')
    print('___________Mean___________')
    extract_pValue(dicMean)

    all[item] = dicMean
    allNone[item] = dicMeanNone
    allCached[item] = dicMeanCached
    allLog[item] = dicMeanLog



# df_new = df_new.iloc[:, 1:]
df = df_new
df.insert(0,'Unnamed: 0',new_cols)
# df.iloc[:, 0] = new_cols
df.set_index(df.columns[0], inplace=True)
df = df.applymap(lambda x: x.replace("'", '"') if isinstance(x, str) else x)

# Convert the strings to dictionaries
df = df.applymap(lambda x: json.loads(x) if isinstance(x, str) else x if pd.notnull(x) else np.nan)
medianprops = dict(linestyle='-', linewidth=1.5, color='black')

for item in df.columns[2:3]:
    dicMean = {}
    dicMedian = {}
    dicStd = {}
    # fig, ax = plt.subplots()
    a = df[item]
    for index, row in a.items():
        for stat in ['mean', 'median', 'std']:
            if stat == 'mean':
                if index not in dicMean:
                    dicMean[index] = [row[stat]]
                else:
                    dicMean[index].append(row[stat])

    dicMeanCached = dicMean.copy()
    for key in dicMean:
        data = dicMean[key]
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        data = np.array(data)
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        mask = (data >= lower_bound) & (data <= upper_bound)
        dicMean[key] = data[mask].tolist()
        # dicMean[key] = np.mean(data[mask].tolist())



    x_values_set = natsort.natsorted(list(dicMean.keys()))
    # x_values_set = ['HC', 'Premanifest', 'Mild', 'Severe']
    float_list = [float(s) for s in x_values_set]
    indexes = np.asarray(float_list)
    indexes[0] = 40
    indexes[-1] = 160

    x_values_set_ticks = x_values_set.copy()
    x_values_set_ticks[0] = 'HC'
    x_values_set_ticks[-1] = 'HGPS'

    dataMean, dataMeanNone, dataMeanCached, dataMeanLog, dataMedian, dataStd = [], [], [], [], [], []
    for x_value in x_values_set:
        currValue = dicMean[x_value]
        dataMeanNone.append(currValue)
        # currValue = np.log(currValue)
        dataMeanLog.append(np.log(dicMean[x_value]))
        # res = [((x - mean[item]) / std[item]) for x in currValue]
        res = [(x / mean[item]) for x in currValue]
        dataMeanCached.append(dicMeanCached[x_value])
        dataMean.append(res)


df_migration = pd.read_excel(r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\MigrationData.xlsx')
Dp_feature_all = df_migration.Dp
Dtot_feature_all = df_migration.Dtot
Pp_feature_all = df_migration.Pp
Psi_feature_all = df_migration.Psi
Pnp_feature_all = df_migration.Pnp
Dnp_feature_all = df_migration.Dnp
MSD_feature_all = df_migration.MSD
Sp_feature_all = df_migration.Sp
Snp_feature_all = df_migration.Snp

fpr_boot = dict()
tpr_boot = dict()
interp_tpr_boot = dict()
roc_auc_boot = dict()
total_bootstrapped = 10000
for boot in range(total_bootstrapped):
    Dp_CAP_all = np.round(df_migration.CAP)
    Dp_CAP = np.unique(Dp_CAP_all)
    Dp_feature_Train = []
    Dp_feature_Test = []
    Dtot_feature_Train = []
    Dtot_feature_Test = []
    Pp_feature_Train = []
    Pp_feature_Test = []
    Psi_feature_Train = []
    Psi_feature_Test = []
    Pnp_feature_Train = []
    Pnp_feature_Test = []
    Dnp_feature_Train = []
    Dnp_feature_Test = []
    MSD_feature_Train = []
    MSD_feature_Test = []
    Sp_feature_Train = []
    Sp_feature_Test = []
    Snp_feature_Train = []
    Snp_feature_Test = []
    split_percentage = 0.6

    for cap_value in Dp_CAP:
        index_map = np.where(Dp_CAP_all == cap_value)[0]
        np.random.shuffle(index_map)
        split_index = int(len(index_map) * split_percentage)
        train_indices = index_map[:split_index]
        test_indices = index_map[split_index:]

        train_mean = np.mean(Dp_feature_all[train_indices])
        train_std = np.std(Dp_feature_all[train_indices])
        train_median = np.median(Dp_feature_all[train_indices])
        test_mean = np.mean(Dp_feature_all[test_indices])
        test_std = np.std(Dp_feature_all[test_indices])
        test_median = np.median(Dp_feature_all[test_indices])
        Dp_feature_Train.append([train_mean, train_std, train_median])
        Dp_feature_Test.append([test_mean, test_std, test_median])

        n_samples = math.ceil(0.3*len(test_indices))
        test_indices_bootstrap = np.random.choice(test_indices, n_samples, replace=False)
        train_mean = np.mean(Dtot_feature_all[train_indices])
        train_std = np.std(Dtot_feature_all[train_indices])
        train_median = np.median(Dtot_feature_all[train_indices])
        test_mean = np.mean(Dtot_feature_all[test_indices_bootstrap])
        test_std = np.std(Dtot_feature_all[test_indices_bootstrap])
        test_median = np.median(Dtot_feature_all[test_indices_bootstrap])
        Dtot_feature_Train.append([train_mean, train_std, train_median])
        Dtot_feature_Test.append([test_mean, test_std, test_median])

        train_mean = np.mean(Pp_feature_all[train_indices])
        train_std = np.std(Pp_feature_all[train_indices])
        train_median = np.median(Pp_feature_all[train_indices])
        test_mean = np.mean(Pp_feature_all[test_indices])
        test_std = np.std(Pp_feature_all[test_indices])
        test_median = np.median(Pp_feature_all[test_indices])
        Pp_feature_Train.append([train_mean, train_std, train_median])
        Pp_feature_Test.append([test_mean, test_std, test_median])

        train_mean = np.mean(Psi_feature_all[train_indices])
        train_std = np.std(Psi_feature_all[train_indices])
        train_median = np.median(Psi_feature_all[train_indices])
        test_mean = np.mean(Psi_feature_all[test_indices])
        test_std = np.std(Psi_feature_all[test_indices])
        test_median = np.median(Psi_feature_all[test_indices])
        Psi_feature_Train.append([train_mean, train_std, train_median])
        Psi_feature_Test.append([test_mean, test_std, test_median])

        train_mean = np.mean(Pnp_feature_all[train_indices])
        train_std = np.std(Pnp_feature_all[train_indices])
        train_median = np.median(Pnp_feature_all[train_indices])
        test_mean = np.mean(Pnp_feature_all[test_indices])
        test_std = np.std(Pnp_feature_all[test_indices])
        test_median = np.median(Pnp_feature_all[test_indices])
        Pnp_feature_Train.append([train_mean, train_std, train_median])
        Pnp_feature_Test.append([test_mean, test_std, test_median])

        train_mean = np.mean(Dnp_feature_all[train_indices])
        train_std = np.std(Dnp_feature_all[train_indices])
        train_median = np.median(Dnp_feature_all[train_indices])
        test_mean = np.mean(Dnp_feature_all[test_indices])
        test_std = np.std(Dnp_feature_all[test_indices])
        test_median = np.median(Dnp_feature_all[test_indices])
        Dnp_feature_Train.append([train_mean, train_std, train_median])
        Dnp_feature_Test.append([test_mean, test_std, test_median])

        train_mean = np.mean(MSD_feature_all[train_indices])
        train_std = np.std(MSD_feature_all[train_indices])
        train_median = np.median(MSD_feature_all[train_indices])
        test_mean = np.mean(MSD_feature_all[test_indices])
        test_std = np.std(MSD_feature_all[test_indices])
        test_median = np.median(MSD_feature_all[test_indices])
        MSD_feature_Train.append([train_mean, train_std, train_median])
        MSD_feature_Test.append([test_mean, test_std, test_median])

        train_mean = np.mean(Sp_feature_all[train_indices])
        train_std = np.std(Sp_feature_all[train_indices])
        train_median = np.median(Sp_feature_all[train_indices])
        test_mean = np.mean(Sp_feature_all[test_indices])
        test_std = np.std(Sp_feature_all[test_indices])
        test_median = np.median(Sp_feature_all[test_indices])
        Sp_feature_Train.append([train_mean, train_std, train_median])
        Sp_feature_Test.append([test_mean, test_std, test_median])

        train_mean = np.mean(Snp_feature_all[train_indices])
        train_std = np.std(Snp_feature_all[train_indices])
        train_median = np.median(Snp_feature_all[train_indices])
        test_mean = np.mean(Snp_feature_all[test_indices])
        test_std = np.std(Snp_feature_all[test_indices])
        test_median = np.median(Snp_feature_all[test_indices])
        Snp_feature_Train.append([train_mean, train_std, train_median])
        Snp_feature_Test.append([test_mean, test_std, test_median])
    Dp_CAP[-1] = 201.

    ch4_feature_total_area_sum = value_lists[10]
    total_area_sum_feature_train = []
    total_area_sum_feature_test = []
    ch1_feature_ratio_width_to_length_mean = value_lists[5]
    ratio_width_to_length_mean_feature_train = []
    ratio_width_to_length_mean_feature_test = []
    biological_CAP = np.round(grouped.index.to_list())

    label_dict = defaultdict(list)
    for label, lst in zip(biological_CAP, ch4_feature_total_area_sum):
        label_dict[label].extend(lst)
    ch4_feature_total_area_sum = list(label_dict.values())
    label_dict = defaultdict(list)
    for label, lst in zip(biological_CAP, ch1_feature_ratio_width_to_length_mean):
        label_dict[label].extend(lst)
    ch1_feature_ratio_width_to_length_mean = list(label_dict.values())

    biological_CAP[-1] = 201.
    biological_CAP = np.unique(biological_CAP)
    for i, cap_value in enumerate(biological_CAP):
        # Gather all instances for current CAP value
        total_area_sum_values = ch4_feature_total_area_sum[i]
        ratio_width_to_length_mean_values = ch1_feature_ratio_width_to_length_mean[i]

        # Calculate the index to split the data
        split_index1 = int(len(total_area_sum_values) * split_percentage)
        split_index2 = int(len(ratio_width_to_length_mean_values) * split_percentage)

        # Split and append values to train and test lists
        train_mean_total_area = np.mean(total_area_sum_values[:split_index1])
        train_std_total_area = np.std(total_area_sum_values[:split_index1])
        train_median_total_area = np.median(total_area_sum_values[:split_index1])
        currentData = total_area_sum_values[split_index1:]
        n_samples = math.ceil(0.3*len(currentData))
        test_indices_bootstrap = np.random.choice(currentData, n_samples, replace=False)
        test_mean_total_area = np.mean(test_indices_bootstrap)
        test_std_total_area = np.std(test_indices_bootstrap)
        test_median_total_area = np.median(test_indices_bootstrap)
        total_area_sum_feature_train.append([train_mean_total_area, train_std_total_area, train_median_total_area])
        total_area_sum_feature_test.append([test_mean_total_area, test_std_total_area, test_median_total_area])

        train_mean_ratio = np.mean(ratio_width_to_length_mean_values[:split_index2])
        train_std_ratio = np.std(ratio_width_to_length_mean_values[:split_index2])
        train_median_ratio = np.median(ratio_width_to_length_mean_values[:split_index2])
        currentData = ratio_width_to_length_mean_values[split_index2:]
        n_samples = math.ceil(0.3*len(currentData))
        test_indices_bootstrap = np.random.choice(currentData, n_samples, replace=False)
        test_mean_ratio = np.mean(test_indices_bootstrap)
        test_std_ratio = np.std(test_indices_bootstrap)
        test_median_ratio = np.median(test_indices_bootstrap)
        ratio_width_to_length_mean_feature_train.append([train_mean_ratio, train_std_ratio, train_median_ratio])
        ratio_width_to_length_mean_feature_test.append([test_mean_ratio, test_std_ratio, test_median_ratio])

    # total_area_sum_feature_train_fixed = []
    # total_area_sum_feature_test_fixed = []
    # ratio_width_to_length_mean_feature_train_fixed = []
    # ratio_width_to_length_mean_feature_test_fixed = []
    # prev = -1
    # use_i = 0
    # for i, cap_value in enumerate(biological_CAP):
    #     if prev == cap_value:
    #         use_i = use_i - 1
    #         total_area_sum_feature_train_fixed[use_i] += total_area_sum_feature_train[i]
    #         total_area_sum_feature_test_fixed[use_i] += total_area_sum_feature_test[i]
    #         ratio_width_to_length_mean_feature_train_fixed[use_i] += ratio_width_to_length_mean_feature_train[i]
    #         ratio_width_to_length_mean_feature_test_fixed[use_i] += ratio_width_to_length_mean_feature_test[i]
    #         total_area_sum_feature_train_fixed[use_i] /= 2
    #         total_area_sum_feature_test_fixed[use_i] /= 2
    #         ratio_width_to_length_mean_feature_train_fixed[use_i] /= 2
    #         ratio_width_to_length_mean_feature_test_fixed[use_i] /= 2
    #     else:
    #         total_area_sum_feature_train_fixed.append(total_area_sum_feature_train[i])
    #         total_area_sum_feature_test_fixed.append(total_area_sum_feature_test[i])
    #         ratio_width_to_length_mean_feature_train_fixed.append(ratio_width_to_length_mean_feature_train[i])
    #         ratio_width_to_length_mean_feature_test_fixed.append(ratio_width_to_length_mean_feature_test[i])
    #     prev = cap_value
    #     use_i = use_i + 1


    form_factor_feature_all = dataMean
    form_factor_CAP = np.array(list(map(float, x_values_set)))
    form_factor_feature_train = []
    form_factor_feature_test = []

    for i, cap_value in enumerate(form_factor_CAP):
        # Gather all instances for current CAP value
        form_factor_values = form_factor_feature_all[i]

        # Calculate the index to split the data
        split_index = int(len(form_factor_values) * split_percentage)

        # Split and append values to train and test lists
        train_mean = np.mean(form_factor_values[:split_index])
        train_std = np.std(form_factor_values[:split_index])
        train_median = np.median(form_factor_values[:split_index])
        currentData = form_factor_values[split_index:]
        n_samples = math.ceil(0.3*len(currentData))
        test_indices_bootstrap = np.random.choice(currentData, n_samples, replace=False)
        test_mean = np.mean(test_indices_bootstrap)
        test_std = np.std(test_indices_bootstrap)
        test_median = np.median(test_indices_bootstrap)
        form_factor_feature_train.append([train_mean, train_std, train_median])
        form_factor_feature_test.append([test_mean, test_std, test_median])
    form_factor_CAP[-1] = 201.

    remove_list = np.array([114., 119., 131., 132., 139.])
    remove_indices = [i for i, cap in enumerate(Dp_CAP) if cap in remove_list]
    Dp_feature_Train = [feat for i, feat in enumerate(Dp_feature_Train) if i not in remove_indices]
    Dp_feature_Test = [feat for i, feat in enumerate(Dp_feature_Test) if i not in remove_indices]
    Dtot_feature_Train = [feat for i, feat in enumerate(Dtot_feature_Train) if i not in remove_indices]
    Dtot_feature_Test = [feat for i, feat in enumerate(Dtot_feature_Test) if i not in remove_indices]
    Pp_feature_Train = [feat for i, feat in enumerate(Pp_feature_Train) if i not in remove_indices]
    Pp_feature_Test = [feat for i, feat in enumerate(Pp_feature_Test) if i not in remove_indices]
    Psi_feature_Train = [feat for i, feat in enumerate(Psi_feature_Train) if i not in remove_indices]
    Psi_feature_Test = [feat for i, feat in enumerate(Psi_feature_Test) if i not in remove_indices]
    Pnp_feature_Train = [feat for i, feat in enumerate(Pnp_feature_Train) if i not in remove_indices]
    Pnp_feature_Test = [feat for i, feat in enumerate(Pnp_feature_Test) if i not in remove_indices]
    Dnp_feature_Train = [feat for i, feat in enumerate(Dnp_feature_Train) if i not in remove_indices]
    Dnp_feature_Test = [feat for i, feat in enumerate(Dnp_feature_Test) if i not in remove_indices]
    MSD_feature_Train = [feat for i, feat in enumerate(MSD_feature_Train) if i not in remove_indices]
    MSD_feature_Test = [feat for i, feat in enumerate(MSD_feature_Test) if i not in remove_indices]
    Sp_feature_Train = [feat for i, feat in enumerate(Sp_feature_Train) if i not in remove_indices]
    Sp_feature_Test = [feat for i, feat in enumerate(Sp_feature_Test) if i not in remove_indices]
    Snp_feature_Train = [feat for i, feat in enumerate(Snp_feature_Train) if i not in remove_indices]
    Snp_feature_Test = [feat for i, feat in enumerate(Snp_feature_Test) if i not in remove_indices]

    remove_indices = [i for i, cap in enumerate(biological_CAP) if cap in remove_list]
    ratio_width_to_length_mean_feature_train = [feat for i, feat in enumerate(ratio_width_to_length_mean_feature_train) if i not in remove_indices]
    ratio_width_to_length_mean_feature_test = [feat for i, feat in enumerate(ratio_width_to_length_mean_feature_test) if i not in remove_indices]
    total_area_sum_feature_train = [feat for i, feat in enumerate(total_area_sum_feature_train) if i not in remove_indices]
    total_area_sum_feature_test = [feat for i, feat in enumerate(total_area_sum_feature_test) if i not in remove_indices]
    remove_indices = [i for i, cap in enumerate(form_factor_CAP) if cap in remove_list]
    form_factor_feature_train = [feat for i, feat in enumerate(form_factor_feature_train) if i not in remove_indices]
    form_factor_feature_test = [feat for i, feat in enumerate(form_factor_feature_test) if i not in remove_indices]

    for i, cap_value in enumerate(remove_list):
        Dp_CAP = np.delete(Dp_CAP, np.where(Dp_CAP == cap_value))
        biological_CAP = np.delete(biological_CAP, np.where(biological_CAP == cap_value))
        form_factor_CAP = np.delete(form_factor_CAP, np.where(form_factor_CAP == cap_value))

    # Split your training and test data
    # train_data = list(zip(Dp_feature_Train, ratio_width_to_length_mean_feature_train,
    #                       total_area_sum_feature_train, form_factor_feature_train))
    # test_data = list(zip(Dp_feature_Test, ratio_width_to_length_mean_feature_test,
    #                      total_area_sum_feature_test, form_factor_feature_test))
    # Concatenate them along the second dimension
    train_data = np.concatenate(
        [Dp_feature_Train, ratio_width_to_length_mean_feature_train, total_area_sum_feature_train, form_factor_feature_train],
        axis=1
    )
    # Concatenate them along the second dimension
    test_data = np.concatenate(
        [Dp_feature_Test, ratio_width_to_length_mean_feature_test, total_area_sum_feature_test, form_factor_feature_test],
        axis=1
    )
    # Rescale the data using MinMaxScaler
    scaler = MinMaxScaler()
    train_data_scaled = scaler.fit_transform(train_data)
    test_data_scaled = scaler.transform(test_data)  # apply same scaling to test data

    # Classification Task
    Dp_CAP_HC = Dp_CAP == 0
    Dp_CAP_Pre = (Dp_CAP > 1) & (Dp_CAP < 90)
    Dp_CAP_Mild = (Dp_CAP >= 90) & (Dp_CAP < 112)
    Dp_CAP_Severe = (Dp_CAP >= 112) & (Dp_CAP < 150)
    Dp_CAP_HGPS = Dp_CAP == 201
    Dp_CAP_Discrete = np.full(Dp_CAP.shape, '', dtype='<U10')
    Dp_CAP_Discrete[Dp_CAP_HC] = 'HC'
    Dp_CAP_Discrete[Dp_CAP_Pre] = 'Premanifest'
    Dp_CAP_Discrete[Dp_CAP_Mild] = 'Mild'
    Dp_CAP_Discrete[Dp_CAP_Severe] = 'Severe'
    Dp_CAP_Discrete[Dp_CAP_HGPS] = 'HGPS'



    def bootstrap_iteration(data, labels, model, n_classes):
        # resample with replacement
        indices = np.arange(len(data))
        resample_indices = resample(indices)
        data_resampled = data[resample_indices]
        labels_resampled = labels[resample_indices]

        # score model on resampled data
        y_score = model.predict(data_resampled)

        # compute ROC curve for each class
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(labels_resampled[:, i], y_score[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

        return fpr, tpr, roc_auc

    path = r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\Models_CLassifier'
    os.makedirs(path, exist_ok=True)

    features_train = [Dp_feature_Train, ratio_width_to_length_mean_feature_train,
                      total_area_sum_feature_train, form_factor_feature_train,
                      Dtot_feature_Train, Dnp_feature_Train, MSD_feature_Train]

    features_test = [Dp_feature_Test, ratio_width_to_length_mean_feature_test,
                     total_area_sum_feature_test, form_factor_feature_test,
                     Dtot_feature_Test, Dnp_feature_Test, MSD_feature_Test]

    # Get all combinations of features
    combs = chain(*map(lambda x: combinations(range(len(features_train)), x), range(1, len(features_train) + 1)))

    # for comb in combs:
    comb = (1, 2, 3, 4)
    train_data = np.concatenate([features_train[i] for i in comb], axis=1)
    test_data = np.concatenate([features_test[i] for i in comb], axis=1)

    scaler = MinMaxScaler()
    train_data_scaled = scaler.fit_transform(train_data)
    test_data_scaled = scaler.transform(test_data)

    le = LabelEncoder()
    labels = le.fit_transform(Dp_CAP_Discrete)
    model = load_model(
        fr'C:\Users\avivs\PycharmProjects\CellInsights\Classifier_Pipeline\Models_Classifier\model_(1, 2, 3, 4).h5')
    test_labels_bin = label_binarize(labels, classes=np.unique(labels))
    n_classes = test_labels_bin.shape[1]

    # Ensure your model's output is probability distributions
    y_score = model.predict(test_data_scaled)

    test_labels_bin_HD = test_labels_bin[:, 2:]
    y_score_HD = y_score[:, 2:]

    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    common_fpr = np.linspace(0, 1, 100)

    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(test_labels_bin[:, i], y_score[:, i])
        interp_tpr = np.interp(common_fpr, fpr[i], tpr[i])
        # Store the interpolated TPR
        if i not in interp_tpr_boot.keys():
            interp_tpr_boot[i] = [interp_tpr]  # store as list of lists
        else:
            interp_tpr_boot[i].append(interp_tpr)  # append to existing list of lists
        # if i not in fpr_boot.keys():
        #     fpr_boot[i] = fpr[i]
        # else:
        #     try:
        #         fpr_boot[i] += fpr[i]
        #     except ValueError:
        #         print(boot, i)
        # if i not in tpr_boot.keys():
        #     try:
        #         tpr_boot[i] = tpr[i]
        #     except ValueError:
        #         print(boot, i)
        # else:
        #     tpr_boot[i] += tpr[i]
        roc_auc[i] = auc(fpr[i], tpr[i])
        if i not in roc_auc_boot.keys():
            roc_auc_boot[i] = roc_auc[i]
        else:
            roc_auc_boot[i] += roc_auc[i]

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(test_labels_bin.ravel(), y_score.ravel())
    interp_tpr_micro = np.interp(common_fpr, fpr["micro"], tpr["micro"])
    if "micro" not in interp_tpr_boot.keys():
        interp_tpr_boot["micro"] = [interp_tpr_micro]  # store as list of lists
    else:
        interp_tpr_boot["micro"].append(interp_tpr_micro)  # append to existing list of lists
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    if "micro" not in roc_auc_boot.keys():
        roc_auc_boot["micro"] = roc_auc["micro"]
    else:
        roc_auc_boot["micro"] += roc_auc["micro"]

    fpr["microHD"], tpr["microHD"], _ = roc_curve(test_labels_bin_HD.ravel(), y_score_HD.ravel())
    interp_tpr_micro = np.interp(common_fpr, fpr["microHD"], tpr["microHD"])
    if "microHD" not in interp_tpr_boot.keys():
        interp_tpr_boot["microHD"] = [interp_tpr_micro]  # store as list of lists
    else:
        interp_tpr_boot["microHD"].append(interp_tpr_micro)  # append to existing list of lists
    roc_auc["microHD"] = auc(fpr["microHD"], tpr["microHD"])
    if "microHD" not in roc_auc_boot.keys():
        roc_auc_boot["microHD"] = roc_auc["microHD"]
    else:
        roc_auc_boot["microHD"] += roc_auc["microHD"]

# fpr_boot = {k: v / total_bootstrapped for k, v in fpr_boot.items()}
# tpr_boot = {k: v / total_bootstrapped for k, v in tpr_boot.items()}
interp_tpr_boot = {k: np.mean(v, axis=0) for k, v in interp_tpr_boot.items()}
roc_auc_boot = {k: v / total_bootstrapped for k, v in roc_auc_boot.items()}

# Plot all ROC curves
plt.figure()
# plt.plot(fpr_boot["micro"], tpr_boot["micro"],
#          label='micro-average ROC curve (area = {0:0.2f})'
#                ''.format(roc_auc_boot["micro"]),
#          color='deeppink', linestyle=':', linewidth=4)

colors = [(66 / 255, 133 / 255, 244 / 255), (51 / 255, 0, 114 / 255),
            (244 / 255, 180 / 255, 0), (15 / 255, 157 / 255, 88 / 255), (219 / 255, 68 / 255, 55 / 255)]
# for i, color in zip(range(n_classes), colors):
#     plt.plot(fpr_boot[i], tpr_boot[i], color=color, lw=2,
#              label='ROC curve of class {0} (area = {1:0.2f})'
#                    ''.format(le.inverse_transform([i])[0], roc_auc_boot[i]))
plt.plot(common_fpr, interp_tpr_boot["micro"], label='micro-average ROC curve (area = {0:0.2f})'
         ''.format(roc_auc_boot["micro"]), color='deeppink', linestyle=':', linewidth=4)
plt.plot(common_fpr, interp_tpr_boot["microHD"], label='HD-only-micro-average ROC curve (area = {0:0.2f})'
         ''.format(roc_auc_boot["microHD"]), color='yellow', linestyle=':', linewidth=4)

for i, color in zip(range(n_classes), colors):
    plt.plot(common_fpr, interp_tpr_boot[i], color=color, lw=2,
             label='ROC curve of class {0} (area = {1:0.2f})'
                   ''.format(le.inverse_transform([i])[0], roc_auc_boot[i]))

plt.plot([0, 1], [0, 1], 'k--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Multi-Class ROC')
plt.legend(loc="lower right")
savePlots = r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\Inference'
os.makedirs(savePlots, exist_ok=True)
plt.savefig(savePlots + fr'\ROC_Bootstrapped10k.png', bbox_inches='tight', dpi=900)

# plt.show()

    # the number of bootstrap iterations
    # n_iterations = 50
    #
    # # store all fpr, tpr, and auc values from each bootstrap iteration
    # fpr_values = []
    # tpr_values = []
    # auc_values = []
    #
    # for _ in range(n_iterations):
    #     fpr, tpr, roc_auc = bootstrap_iteration(test_data_scaled, test_labels_bin, model, n_classes)
    #     fpr_values.append(fpr)
    #     tpr_values.append(tpr)
    #     auc_values.append(roc_auc)
    #
    # # Define a grid of FPR values (from 0 to 1) to interpolate TPR values.
    # mean_fpr = np.linspace(0, 1, 100)
    #
    # # Empty lists to hold the interpolated TPR values and AUC values for each bootstrap iteration.
    # mean_tpr = []
    # roc_auc = []
    #
    # for i in range(n_iterations):
    #     # For each class
    #     for j in range(n_classes):
    #         # Interpolate the TPR for current bootstrap iteration ROC onto the common FPR grid
    #         interp_tpr = interp(mean_fpr, fpr_values[i][j], tpr_values[i][j])
    #         # Append the interpolated TPR to the list
    #         mean_tpr.append(interp_tpr)
    #         # Calculate and append the AUC for current bootstrap iteration
    #         roc_auc.append(auc(mean_fpr, interp_tpr))
    #
    # # Calculate the mean interpolated TPR across all bootstrap iterations for each point on the FPR grid.
    # mean_tpr = np.mean(mean_tpr, axis=0)
    # # Calculate the mean AUC across all bootstrap iterations.
    # mean_roc_auc = auc(mean_fpr, mean_tpr)
    #
    # # Plot the mean ROC curve.
    # plt.figure()
    # plt.plot(mean_fpr, mean_tpr, label='Mean ROC curve (area = {:.2f})'.format(mean_roc_auc))
    # plt.plot([0, 1], [0, 1], 'k--')
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    # plt.xlabel('False Positive Rate')
    # plt.ylabel('True Positive Rate')
    # plt.title('Bootstrap Mean ROC')
    # plt.legend(loc='lower right')
    # plt.show()
