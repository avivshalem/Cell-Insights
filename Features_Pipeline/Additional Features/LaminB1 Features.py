import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score
import scipy.stats as stats
from scipy.stats import shapiro, kstest, anderson, norm
from statsmodels.graphics.gofplots import qqplot
from sklearn.ensemble import IsolationForest
import os

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
    plt.savefig(savePlots + fr'\{feature}_QQ_Plot.png', bbox_inches='tight')

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


def extract_pValue_return(dic):
    normal = check_normal_dist(dic)
    p_value_dict = {}
    for key1 in dic.keys():
        for key2 in dic.keys():
            if key1 != key2:
                pair_key = f"{key1}_vs_{key2}"
                t_statistic, p_value = stats.ttest_ind(dic[key1], dic[key2])
                u_statistic, p_valueMW = stats.mannwhitneyu(dic[key1], dic[key2])

                if normal:
                    p_value_dict[pair_key] = p_value
                else:
                    p_value_dict[pair_key] = p_valueMW
                print(pair_key, p_value_dict[pair_key])
    return p_value_dict

# Load the Excel file into a Pandas DataFrame
file = 'PlateResults - LaminB+a actinin- S exp4 05.07.22__2022-07-12T10_09_44-Measurement 1'
df1 = pd.read_excel(rf'\LaminB_modified\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- S exp4 05.07.22 (9-10)__2022-07-12T11_44_17-Measurement 1'
df2 = pd.read_excel(rf'\LaminB_modified\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- R exp4 05.07.22__2022-07-17T09_42_44-Measurement 1'
df3 = pd.read_excel(rf'\LaminB_modified\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- exp8 12.02.23__2023-02-20T14_52_23-Measurement 1'
df4 = pd.read_excel(rf'\LaminB_modified\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- exp7 06.02.23__2023-02-19T14_31_42-Measurement 1'
df5 = pd.read_excel(rf'\LaminB_modified\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- exp6  05.02.23__2023-02-19T12_14_21-Measurement 1'
df6 = pd.read_excel(rf'\LaminB_modified\{file}.xlsx')
file = 'PlateResults - LaminB+a actinin- exp5 plate1 08.01.23__2023-01-12T15_50_14-Measurement 1'
df7 = pd.read_excel(rf'\LaminB_modified\{file}.xlsx')
file = 'HGPS PlateResults'
df8 = pd.read_excel(rf'\LaminB_modified\{file}.xlsx')
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
f_names = []
for fi, feature in enumerate(features):
    # if "ir_chanel_4_total Area" not in feature and "chanel_1 Ratio Width to Length" not in feature:
    # if "chanel_1 Ratio Width to Length" not in feature:
    if "chanel_4" not in feature:
        continue
    else:
        f_names.append(feature)
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


savePlots = rf'\LaminB1'
os.makedirs(savePlots, exist_ok=True)
# Group the data by the values in column A and aggregate the other columns by taking their mean
medianprops = dict(linestyle='-', linewidth=1.5, color='black')

features = df.columns[4:-6].to_list()
# counter = 0
for fi, feature in enumerate(features):
    # if "ir_chanel_4_total Area" not in feature and "chanel_1 Ratio Width to Length" not in feature:
    # if "chanel_1 Ratio Width to Length" not in feature:
    if "chanel_4" not in feature:
        continue
    changer = False
    if "ir_chanel" in feature:
        changer = True
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

    for key, value in list(dicRaw.items()):
        # Check if the value is empty or None
        if len(value) == 0:
            # Remove the key/value pair from the dictionary
            del dicRaw[key]
    for key, value in list(dicOutliers.items()):
        # Check if the value is empty or None
        if len(value) == 0:
            # Remove the key/value pair from the dictionary
            del dicOutliers[key]
    for key, value in list(dic.items()):
        # Check if the value is empty or None
        if len(value) == 0:
            # Remove the key/value pair from the dictionary
            del dic[key]
    p_value_dict = extract_pValue_return(dic)
    # extract_pValue(dic)

    indexesOG = grouped.index.to_list()
    indexesOG = np.round(indexesOG)
    indexes = grouped.index.to_list()
    indexes[0] = 40
    indexes[-1] = 161
    indexes = np.round(indexes)

    # fig11C, axsC = plt.subplots(int(len(groupedRaw.values) / 4), 4, figsize=(6, 8))
    # plt.subplots_adjust(hspace=0.5)
    # fig11C.suptitle('Distribution of Data Raw')
    # axsC = axsC.flatten()
    # for i, lst in enumerate(groupedRaw.values):
    #     axsC[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axsC[i].set_title(indexes[i])
    #     axsC[i].set_xlim(-1.2*Xlim, 1.2*Xlim)
    #     # axs2[i].set_ylim(0, 40)
    # axsC[-1].set_xlabel('Value')
    # fig11C.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    #
    # fig11N, axsN = plt.subplots(int(len(groupedOutliers.values) / 4), 4, figsize=(6, 8))
    # plt.subplots_adjust(hspace=0.5)
    # fig11N.suptitle('Distribution of Data Outlier Removed')
    # axsN = axsN.flatten()
    # for i, lst in enumerate(groupedOutliers.values):
    #     axsN[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axsN[i].set_title(indexes[i])
    #     axsN[i].set_xlim(-1.2*Xlim, 1.2*Xlim)
    #     # axs2[i].set_ylim(0, 40)
    # axsN[-1].set_xlabel('Value')
    # fig11N.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    #
    # fig11L, axsL = plt.subplots(int(len(groupedLog.values) / 4), 4, figsize=(6, 8))
    # plt.subplots_adjust(hspace=0.5)
    # fig11L.suptitle('Distribution of Log of Data')
    # axsL = axsL.flatten()
    # for i, lst in enumerate(groupedLog.values):
    #     axsL[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axsL[i].set_title(indexes[i])
    #     axsL[i].set_xlim(-1.2*Xlim, 1.2*Xlim)
    #     # axs2[i].set_ylim(0, 40)
    # axsL[-1].set_xlabel('Value')
    # fig11L.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    #
    # fig11, axs = plt.subplots(int(len(grouped.values) / 4), 4, figsize=(6, 8))
    # fig11.suptitle('Distribution of Data Mean Normalization')
    # axs = axs.flatten()
    # plt.subplots_adjust(hspace=0.5)
    # for i, lst in enumerate(grouped.values):
    #     axs[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axs[i].set_title(indexes[i])
    #     axs[i].set_xlim(-1.2*Xlim, 1.2*Xlim)
    #     # axs[i].set_ylim(0, 40)
    # axs[-1].set_xlabel('Value')
    # fig11.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(6, 8))

    # valueUse = (grouped.values - mean) / std
    valueUse = (grouped.values) / mean
    bp = ax.boxplot(valueUse, positions=indexes, patch_artist=True, medianprops=medianprops)
    ax.set_xlabel('CAP')
    ax.set_ylabel('Normalized Value')
    ax.set_title(feature)
    # ax.set_xlim([-1000, 160*1000])
    indexes_list = indexesOG.tolist()
    indexes_str = [str(i) for i in indexes_list]
    indexes_str[0] = "HC"
    indexes_str[-1] = "HGPS"

    ax.set_xticklabels(indexes_str, fontsize=5, rotation=70)
    for i, box in enumerate(bp['boxes']):
        if indexes[i] == 40:
            box.set(facecolor='blue')
        elif 40 < indexes[i] < 90:
            box.set(facecolor='green')
        elif 90 <= indexes[i] < 114:
            box.set(facecolor='orange')
        elif indexes[i] == 161:
            box.set(facecolor='purple')
        else:
            box.set(facecolor='red')
    delList= []
    for val in range(len(valueUse)):
        if len(valueUse[val]) == 0:
            delList.append(val)
    count = 0
    for de in delList:
        valueUse = np.delete(valueUse, de-count)
        indexes = np.delete(indexes, de-count)
        count += 1
    x = np.array(indexes)
    y = np.array([np.mean(values) for values in valueUse])
    z = np.polyfit(x, y, 2)
    p = np.poly1d(z)
    # ax.plot(x, p(x), color='magenta')
    y_pred = p(x)
    r2 = r2_score(y, y_pred)

    # Display the R-squared value on the plot
    # ax.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax.transAxes, fontsize=12)

    # xticks = np.arra (grouped.index.min(), grouped.index.max() + 1, 10)  # example with ticks every 10
    # ax.set_xticks(xticks)
    order = list(dic.keys())
    box_plot_dataCached, box_plot_dataNone, box_plot_dataLog, box_plot_data = [], [], [] , []
    indexes = [x+1 for x in range(len(order))]
    for key, value in dicRaw.items():
        box_plot_dataCached.append(value)
    for key, value in dicOutliers.items():
        box_plot_dataNone.append(value)
    for key, value in dicLog.items():
        box_plot_dataLog.append(value)
    for key, value in dic.items():
        box_plot_data.append(value)
    # Xlim = max([max(max(box_plot_dataCached)),max(max(box_plot_dataNone)),max(max(box_plot_data))])

    fig22C, axs2C = plt.subplots(len(box_plot_dataCached), 1, figsize=(6, 8))
    plt.subplots_adjust(hspace=0.5)
    fig22C.suptitle(f'{feature} Distribution of Data Raw')
    gather_lst = []
    for i, lst in enumerate(box_plot_dataCached):
        gather_lst.extend(lst)
    for i, lst in enumerate(box_plot_dataCached):
        n, bins, patches = axs2C[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
        axs2C[i].set_title(order[i])
        axs2C[i].set_xlim([min(gather_lst) - 0.2*min(gather_lst), max(gather_lst) + 0.2*max(gather_lst)])
        axs2C[i].set_ylim([0, max(n)+2])
    axs2C[-1].set_xlabel('Value')
    fig22C.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    fig22C.savefig(savePlots + fr'\{feature}_Raw.png', bbox_inches='tight')

    fig22N, axs2N = plt.subplots(len(box_plot_dataNone), 1, figsize=(6, 8))
    plt.subplots_adjust(hspace=0.5)
    fig22N.suptitle(f'{feature} Distribution of Data Outlier Removed')
    gather_lst = []
    for i, lst in enumerate(box_plot_dataNone):
        gather_lst.extend(lst)
    for i, lst in enumerate(box_plot_dataNone):
        n, bins, patches = axs2N[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
        axs2N[i].set_title(order[i])
        axs2N[i].set_xlim([min(gather_lst) - 0.2*min(gather_lst), max(gather_lst) + 0.2*max(gather_lst)])
        axs2N[i].set_ylim([0, max(n)+2])
    axs2N[-1].set_xlabel('Value')
    fig22N.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    fig22N.savefig(savePlots + fr'\{feature}_OL.png', bbox_inches='tight')

    # fig22L, axs2L = plt.subplots(len(box_plot_dataLog), 1, figsize=(6, 8))
    # plt.subplots_adjust(hspace=0.5)
    # fig22L.suptitle(f'{feature} Distribution of Data Log of data')
    # for i, lst in enumerate(box_plot_dataLog):
    #     axs2L[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axs2L[i].set_title(order[i])
    #     axs2L[i].set_xlim(-1.2*Xlim, 1.2*Xlim)
    #     # axs2[i].set_ylim(0, 40)
    # axs2L[-1].set_xlabel('Value')
    # fig22L.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    # fig22L.savefig(savePlots + fr'\{feature}_Log.png', bbox_inches='tight')

    fig22, axs2 = plt.subplots(len(box_plot_data), 1, figsize=(6, 8))
    plt.subplots_adjust(hspace=0.5)
    fig22.suptitle(f'{feature} Distribution of Data Mean Normalization')
    gather_lst = []
    for i, lst in enumerate(box_plot_data):
        gather_lst.extend(lst)
    for i, lst in enumerate(box_plot_data):
        n, bins, patches = axs2[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
        axs2[i].set_title(order[i])
        axs2[i].set_xlim([min(gather_lst) - 0.2*min(gather_lst), max(gather_lst) + 0.2*max(gather_lst)])
        axs2[i].set_ylim([0, max(n)+2])
    axs2[-1].set_xlabel('Value')
    fig22.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    fig22.savefig(savePlots + fr'\{feature}_Mean.png', bbox_inches='tight')

    bp2 = ax2.boxplot(box_plot_data, positions=indexes, patch_artist=True, medianprops=medianprops)
    ax2.set_xticklabels(dic.keys())
    current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    offset = 0
    if "chanel_1" in feature:
        current_max_y = 1.08  # np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    #     tit = 'Cells Width to Length Ratio'
    elif "chanel_4 Gabor" in feature:
        incr = 0.025
        current_max_y = 1.05#np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    #     tit = 'LaminB1 Total Area'
    elif "chanel_4 SER Hole" in feature:
        incr = 0.005
        current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    elif "chanel_4 SER Ridge" in feature:
        incr = 0.01
        current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    elif "chanel_4 Spot" in feature:
        incr = 0.025
        current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    elif "chanel_4 SER Valley" in feature:
        incr = 0.035
        current_max_y = 1.25
    elif "chanel_4_Intensity Mean" in feature:
        incr = 0.075
        current_max_y = 1.5
    elif "chanel_4_Intensity Sum" in feature:
        incr = 0.075
        current_max_y = 1.55
    elif "chanel_4_at Mean" in feature:
        incr = 0.025
        current_max_y = 1.1
    elif "chanel_4_at Area" in feature:
        incr = 0.05
        current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    elif "chanel_4_total Area" in feature:
        incr = 0.085
        current_max_y = 1.55
    elif "chanel_4_normalized" in feature:
        incr = 0.035
        current_max_y = 1.25
    ax2.set_ylabel('Normalized Value')
    # ax2.set_title(feature)
    for i, box in enumerate(bp2['boxes']):
        if indexes[i] == 1:
            box.set(facecolor='blue', alpha=0.5)
        elif indexes[i] == 2:
            box.set(facecolor='green', alpha=0.5)
        elif indexes[i] == 3:
            box.set(facecolor='orange', alpha=0.5)
        elif indexes[i] == 5:
            box.set(facecolor='purple', alpha=0.5)
        else:
            box.set(facecolor='red', alpha=0.5)
    x = np.array(indexes)
    y = np.array([np.mean(values) for values in box_plot_data])
    z = np.polyfit(x, y, 2)
    p = np.poly1d(z)
    # ax2.plot(x, p(x), color='magenta')
    y_pred = p(x)
    r2 = r2_score(y, y_pred)

    # Display the R-squared value on the plot
    # ax2.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax2.transAxes, fontsize=12)

    # plt.show()

# for i, feature in enumerate(features):
#     plt.figure()
#     X, Y, E = remove_outliers(grouped.index, grouped[feature], std[feature])
#     plt.scatter(X[X==0], Y[X==0], label=feature, c='blue')
#     plt.scatter(X[(X>0) & (X<90)], Y[(X>0) & (X<90)], label=feature, c='green')
#     plt.scatter(X[(X>=90) & (X<114)], Y[(X>=90) & (X<114)], label=feature, c='orange')
#     plt.scatter(X[X>=114], Y[X>=114], label=feature, c='red')
#     plt.errorbar(X, Y, yerr=E, fmt='none', capsize=5)
#     # Add labels and legend to the plot
#     plt.xlabel('CAP')
#     plt.ylabel('Feature Value')
#     plt.title(feature)
#     # Show the plot

    # Annotate p-values with stars
    for i, key1 in enumerate(dic.keys()):
        for j, key2 in enumerate(dic.keys()):
            if i > j:
                continue
            if key1 != key2:
                pair_key = f"{key1}_vs_{key2}"
                p_value = p_value_dict[pair_key]
                star_annotation = 'N.P.'
                if p_value < 0.001:
                    star_annotation = '***'
                elif p_value < 0.01:
                    star_annotation = '**'
                elif p_value < 0.05:
                    star_annotation = '*'

                if star_annotation:
                    x_line_center = (indexes[i] + indexes[j]) / 2
                    y_line_center = current_max_y + offset

                    # Draw the line
                    ax2.plot([indexes[i], indexes[j]], [y_line_center, y_line_center], color='grey', linestyle='--')

                    # Annotate text
                    ax2.annotate(star_annotation, xy=(x_line_center, y_line_center),
                                 xytext=(0, 2), textcoords='offset points', ha='center', fontsize=6)

                    # Increase the height for the next line

                    offset += current_max_y * incr  # Increment by an additional 5% of the max y-value

    fig.savefig(savePlots + fr'\{feature}_Box.png', bbox_inches='tight', dpi=900)
    # plt.show()

    fig2, ax3 = plt.subplots(1, 1)

    bp2 = ax3.boxplot(box_plot_data, positions=indexes, patch_artist=True, medianprops=medianprops)
    ax3.set_xticklabels(dic.keys())
    current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    if "chanel_1" in feature:
        current_max_y = 1.08  # np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    #     tit = 'Cells Width to Length Ratio'
    elif "chanel_4 Gabor" in feature:
        incr = 0.025
        current_max_y = 1.05#np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    #     tit = 'LaminB1 Total Area'
    elif "chanel_4 SER Hole" in feature:
        incr = 0.01
        current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    elif "chanel_4 SER Ridge" in feature:
        incr = 0.01
        current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    elif "chanel_4 Spot" in feature:
        incr = 0.025
        current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    elif "chanel_4 Valley" in feature:
        incr = 0.035
        current_max_y = 1.25
    elif "chanel_4_Intensity Mean" in feature:
        incr = 0.05
        current_max_y = 1.5
    elif "chanel_4_Intensity Sum" in feature:
        incr = 0.05
        current_max_y = 1.55
    elif "chanel_4_at Mean" in feature:
        incr = 0.025
        current_max_y = 1.1
    elif "chanel_4_at Area" in feature:
        incr = 0.05
        current_max_y = np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    elif "chanel_4_total Area" in feature:
        incr = 0.05
        current_max_y = 1.55
    elif "chanel_4_normalized" in feature:
        incr = 0.035
        current_max_y = 1.25
    offset = 0
    ax3.set_ylabel('Normalized Value')
    # ax2.set_title(feature)
    for i, box in enumerate(bp2['boxes']):
        if indexes[i] == 1:
            box.set(facecolor='blue', alpha=0.5)
        elif indexes[i] == 2:
            box.set(facecolor='green', alpha=0.5)
        elif indexes[i] == 3:
            box.set(facecolor='orange', alpha=0.5)
        elif indexes[i] == 5:
            box.set(facecolor='purple', alpha=0.5)
        else:
            box.set(facecolor='red', alpha=0.5)
    x = np.array(indexes)
    y = np.array([np.median(values) for values in box_plot_data])
    z = np.polyfit(x, y, 2)
    p = np.poly1d(z)
    # ax3.plot(x, p(x), color='magenta')
    y_pred = p(x)
    r2 = r2_score(y, y_pred)

    # ax3.set_title(tit)
    # Display the R-squared value on the plot
    # ax3.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax3.transAxes, fontsize=12)

    # Annotate p-values with stars
    for i, key1 in enumerate(dic.keys()):
        for j, key2 in enumerate(dic.keys()):
            if i > j:
                continue
            if key1 != key2:
                pair_key = f"{key1}_vs_{key2}"
                p_value = p_value_dict[pair_key]
                star_annotation = 'N.P.'
                if p_value < 0.001:
                    star_annotation = '***'
                elif p_value < 0.01:
                    star_annotation = '**'
                elif p_value < 0.05:
                    star_annotation = '*'

                if star_annotation:
                    x_line_center = (indexes[i] + indexes[j]) / 2
                    y_line_center = current_max_y + offset

                    # Draw the line
                    ax3.plot([indexes[i], indexes[j]], [y_line_center, y_line_center], color='grey', linestyle='--')

                    # Annotate text
                    ax3.annotate(star_annotation, xy=(x_line_center, y_line_center),
                                 xytext=(0, 2), textcoords='offset points', ha='center', fontsize=6)

                    # Increase the height for the next line

                    if changer:
                        incr = 0.1
                    offset += current_max_y * incr  # Increment by an additional 5% of the max y-value

    fig2.savefig(savePlots + fr'\{feature}_Box_Grouped.png', bbox_inches='tight', dpi=900)

pass
# print(counter)