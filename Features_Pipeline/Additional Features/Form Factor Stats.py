import os
import pandas as pd
import matplotlib.pyplot as plt
import json
import numpy as np
import scipy.stats as stats
import natsort
from sklearn.metrics import r2_score
from scipy.stats import shapiro, kstest, anderson, norm
from statsmodels.graphics.gofplots import qqplot

from sklearn.ensemble import IsolationForest
import seaborn as sns

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
    plt.savefig(savePlots + fr'\{item}_QQ_Plot.png', bbox_inches='tight')

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
    return p_value_dict

CounterXL = 0
df1 = pd.read_excel(r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\output.xlsx')
df2 = pd.read_excel(r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\output707_708_24H_normed.xlsx')
df3 = pd.read_excel(r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\outputHCR.xlsx')
df = pd.concat([df1, df2, df3])
df_new = pd.DataFrame()
savePlots = r'C:\CellInsights\Form_Factor'
os.makedirs(savePlots, exist_ok=True)
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

# x_axis = ['HC', 'Premanifest', 'Mild', 'Severe', 'HGPS']
# y_axis_circ_mean = [0 for i in range(5)]
# y_axis_enclosing_circle_mean = [0 for i in range(5)]
# y_axis_enclosing_ellipse_mean = [0 for i in range(5)]
# y_axis_ratio_mean = [0 for i in range(5)]
# y_axis_ecc_mean = [0 for i in range(5)]
# y_axis_form_mean = [0 for i in range(5)]
# y_axis_circ_median = [0 for i in range(5)]
# y_axis_enclosing_circle_median = [0 for i in range(5)]
# y_axis_enclosing_ellipse_median = [0 for i in range(5)]
# y_axis_ratio_median = [0 for i in range(5)]
# y_axis_ecc_median = [0 for i in range(5)]
# y_axis_form_median = [0 for i in range(5)]
# y_axis_circ_std = [0 for i in range(5)]
# y_axis_enclosing_circle_std = [0 for i in range(5)]
# y_axis_enclosing_ellipse_std = [0 for i in range(5)]
# y_axis_ratio_std = [0 for i in range(5)]
# y_axis_ecc_std = [0 for i in range(5)]
# y_axis_form_std = [0 for i in range(5)]
#
# df.iloc[:, 0] = new_cols
# df.set_index(df.columns[0], inplace=True)
# for index, row in df.iterrows():
#     for jindex, column in enumerate(df.columns):
#         value = row[column]
#         json_str = value.replace("'", "\"")
#         # Convert JSON string to dictionary
#         dict_obj = json.loads(json_str)
#         if index == 'HC':
#             if column == 'Circularity':
#                 y_axis_circ_mean[0] = dict_obj['mean']
type_cols = []
for col in range(len(new_cols)):
    if 'HC' in new_cols[col]:
        type_cols.append('HC')
    elif 'Premanifest' in new_cols[col]:
        type_cols.append('Premanifest')
    if 'Mild' in new_cols[col]:
        type_cols.append('Mild')
    if 'Severe' in new_cols[col]:
        type_cols.append('Severe')
    if 'HGPS' in new_cols[col]:
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
            # if stat == 'median':
            #     if index not in dicMedian:
            #         dicMedian[index] = [row[stat]]
            #     else:
            #         dicMedian[index].append(row[stat])
            # if stat == 'std':
            #     if index not in dicStd:
            #         dicStd[index] = [row[stat]]
            #     else:
            #         dicStd[index].append(row[stat])
    # outlier removal
    dicMeanCached = dicMean.copy()
    for key in dicMean:
        data = np.asarray(dicMean[key])
        # model = IsolationForest(contamination=0.33)
        # model.fit(data.reshape(-1, 1))
        # outliers = model.predict(data.reshape(-1, 1))
        # data_clean = data[outliers == 1]
        # dicMean[key] = data_clean.tolist()
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        data = np.array(data)
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        mask = (data >= lower_bound) & (data <= upper_bound)
        dicMean[key] = data[mask].tolist()
        # dicMean[key] = np.mean(data[mask].tolist())
    # for key in dicMedian:
    #     data = dicMedian[key]
    #     q1 = np.percentile(data, 25)
    #     q3 = np.percentile(data, 75)
    #     iqr = q3 - q1
    #     data = np.array(data)
    #     lower_bound = q1 - 1.5 * iqr
    #     upper_bound = q3 + 1.5 * iqr
    #     mask = (data >= lower_bound) & (data <= upper_bound)
    #     dicMedian[key] = data[mask].tolist()
    # for key in dicStd:
    #     data = dicStd[key]
    #     q1 = np.percentile(data, 25)
    #     q3 = np.percentile(data, 75)
    #     iqr = q3 - q1
    #     data = np.array(data)
    #     lower_bound = q1 - 1.5 * iqr
    #     upper_bound = q3 + 1.5 * iqr
    #     mask = (data >= lower_bound) & (data <= upper_bound)
    #     dicStd[key] = data[mask].tolist()

    # print('___________Median___________')
    # extract_pValue(dicMedian)
    # print('___________STD___________')
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
        #outers alll
        # data = np.asarray(dicMean[key])
        # model = IsolationForest(contamination='auto')
        # model.fit(data.reshape(-1, 1))
        # outliers = model.predict(data.reshape(-1, 1))
        # data_clean = data[outliers == 1]
        # dicMean[key] = data_clean.tolist()
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

for item in df.columns:
    # if item != 'Form Factor':
    #     continue
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
            # if stat == 'median':
            #     if index not in dicMedian:
            #         dicMedian[index] = [row[stat]]
            #     else:
            #         dicMedian[index].append(row[stat])
            # if stat == 'std':
            #     if index not in dicStd:
            #         dicStd[index] = [row[stat]]
            #     else:
            #         dicStd[index].append(row[stat])


    # outlier removal
    # for key in dicMean:
    #     data = np.asarray(dicMean[key])
        # model = IsolationForest(contamination=0.33)
        # model.fit(data.reshape(-1, 1))
        # outliers = model.predict(data.reshape(-1, 1))
        # data_clean = data[outliers == 1]
        # dicMean[key] = data_clean.tolist()
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


    # for key in dicMedian:
    #     data = dicMedian[key]
    #     q1 = np.percentile(data, 25)
    #     q3 = np.percentile(data, 75)
    #     iqr = q3 - q1
    #     data = np.array(data)
    #     lower_bound = q1 - 1.5 * iqr
    #     upper_bound = q3 + 1.5 * iqr
    #     mask = (data >= lower_bound) & (data <= upper_bound)
    #     dicMedian[key] = data[mask].tolist()
    # for key in dicStd:
    #     data = dicStd[key]
    #     q1 = np.percentile(data, 25)
    #     q3 = np.percentile(data, 75)
    #     iqr = q3 - q1
    #     data = np.array(data)
    #     lower_bound = q1 - 1.5 * iqr
    #     upper_bound = q3 + 1.5 * iqr
    #     mask = (data >= lower_bound) & (data <= upper_bound)
    #     dicStd[key] = data[mask].tolist()


    # lower, higher = 25, 75
    # dicMean = alignment_values_to_all(dicMean, 'HC')
    # dicMedian = alignment_values_to_all(dicMedian, 'HC')
    # dicStd = alignment_values_to_all(dicStd, 'HC')
    # dicMean = alignment_values(dicMean, 'Severe', 'HGPS', lower, higher)
    # dicMedian = alignment_values(dicMedian, 'Severe', 'HGPS', lower, higher)
    # dicStd = alignment_values(dicStd, 'Severe', 'HGPS', lower, higher)
    #
    # dicMean = alignment_values(dicMean, 'Mild', 'Severe', lower, higher)
    # dicMedian = alignment_values(dicMedian, 'Mild', 'Severe', lower, higher)
    # dicStd = alignment_values(dicStd, 'Mild', 'Severe', lower, higher)
    #
    # dicMean = alignment_values(dicMean, 'Premanifest', 'Mild', lower, higher)
    # dicMedian = alignment_values(dicMedian, 'Premanifest', 'Mild', lower, higher)
    # dicStd = alignment_values(dicStd, 'Premanifest', 'Mild', lower, higher)

    # dicMean = alignment_values(dicMean, 'HC', 'Premanifest', lower, higher)
    # dicMedian = alignment_values(dicMedian, 'HC', 'Premanifest', lower, higher)
    # dicStd = alignment_values(dicStd, 'HC', 'Premanifest', lower, higher)
    # # #

    # #

    #


    # lower, higher = 85, 25
    # dicMean = alignment_values(dicMean, 'Mild', 'HC', lower, higher)
    # dicMedian = alignment_values(dicMedian, 'Mild', 'HC', lower, higher)
    # dicStd = alignment_values(dicStd, 'Mild', 'HC', lower, higher)
    # lower, higher = 95, 25
    # dicMean = alignment_values(dicMean, 'Severe', 'HC', lower, higher)
    # dicMedian = alignment_values(dicMedian, 'Severe', 'HC', lower, higher)
    # dicStd = alignment_values(dicStd, 'Severe', 'HC', lower, higher)
    # lower, higher = 99, 25
    # dicMean = alignment_values(dicMean, 'HGPS', 'HC', lower, higher)
    # dicMedian = alignment_values(dicMedian, 'HGPS', 'HC', lower, higher)
    # dicStd = alignment_values(dicStd, 'HGPS', 'HC', lower, higher)

    # print(f'***************************{item}************************')
    # x_valuesMean = list(dicMean.keys())
    # print('___________Mean___________')
    # extract_pValue(dicMean)
    # print('___________Median___________')
    # extract_pValue(dicMedian)
    # print('___________STD___________')
    # extract_pValue(dicStd)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))

    ax1.set_title(f"Box Plot Mean for {item}")
    ax1.set_xlabel("CAP")
    ax1.set_ylabel("Feature Value")
    # fig1, (ax11, ax22) = plt.subplots(2, 1, figsize=(6, 8))
    # ax11.set_title(f"Distribution plot for {item} by CAP")
    # ax11.set_xlabel("CAP")
    # ax11.set_ylabel("Occurrences")
    # ax22.set_title(f"Distribution plot for {item}")
    # ax22.set_ylabel("Occurrences")

    # x_valuesMedian = list(dicMedian.keys())
    # ax2.set_title(f"Box Plot Median for {item}")
    # ax2.set_xlabel("X values")
    # ax2.set_ylabel("Y values")
    # x_valuesStd = list(dicStd.keys())
    # ax3.set_title(f"Box Plot STD for {item}")
    # ax3.set_xlabel("X values")
    # ax3.set_ylabel("Y values")

    # Loop through the x_values and plot the corresponding y_values
    # for x_value in x_valuesMean:
    #     ax1.scatter([x_value] * len(dicMean[x_value]), dicMean[x_value])
    # for x_value in x_valuesMedian:
    #     ax2.scatter([x_value] * len(dicMedian[x_value]), dicMedian[x_value])
    # for x_value in x_valuesStd:
    #     ax3.scatter([x_value] * len(dicStd[x_value]), dicStd[x_value])
    # fig.subplots_adjust(hspace=0.5)
    # plt.show()
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
        # ax1.scatter([x_value] * len(dicMean[x_value]), dicMean[x_value])
        # ax1.scatter(x_value, dicMean[x_value])
    # for x_value in x_values_set:
    #     dataMedian.append(dicMedian[x_value])
    #     # ax2.scatter([x_value] * len(dicMedian[x_value]), dicMedian[x_value])
    #     # ax2.scatter(x_value, dicMedian[x_value])
    # for x_value in x_values_set:
    #     dataStd.append(dicStd[x_value])
    #     # ax3.scatter([x_value] * len(dicStd[x_value]), dicStd[x_value])
    #     # ax3.scatter(x_value, dicStd[x_value])
    # for i, group in enumerate(dataMean):
    #     # Combine all the values into a single list
    #     values = [val for val in group]
    #     sns.displot(values, ax=ax11)

    # fig11C, axsC = plt.subplots(int(len(dataMeanCached) / 4), 4, figsize=(6, 8))
    # plt.subplots_adjust(hspace=0.5)
    # fig11C.suptitle(f'{item} Distribution of Data Raw')
    # axsC = axsC.flatten()
    # for i, lst in enumerate(dataMeanCached):
    #     axsC[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axsC[i].set_title(indexes[i])
    #     axsC[i].set_xlim(-2, 2)
    #     # axs2[i].set_ylim(0, 40)
    # axsC[-1].set_xlabel('Value')
    # fig11C.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    #
    # fig11N, axsN = plt.subplots(int(len(dataMeanNone)/4), 4, figsize=(6, 8))
    # plt.subplots_adjust(hspace=0.5)
    # fig11N.suptitle(f'{item} Distribution of Data Outlier Removed')
    # axsN = axsN.flatten()
    # for i, lst in enumerate(dataMeanNone):
    #     axsN[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axsN[i].set_title(indexes[i])
    #     axsN[i].set_xlim(-2, 2)
    #     # axs2[i].set_ylim(0, 40)
    # axsN[-1].set_xlabel('Value')
    # fig11N.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    #
    # fig11L, axsL = plt.subplots(int(len(dataMeanLog) / 4), 4, figsize=(6, 8))
    # plt.subplots_adjust(hspace=0.5)
    # fig11L.suptitle(f'{item} Distribution of Log of Data')
    # axsL = axsL.flatten()
    # for i, lst in enumerate(dataMeanLog):
    #     axsL[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axsL[i].set_title(indexes[i])
    #     axsL[i].set_xlim(-2, 2)
    #     # axs2[i].set_ylim(0, 40)
    # axsL[-1].set_xlabel('Value')
    # fig11L.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    #
    # fig11, axs = plt.subplots(int(len(dataMean)/4), 4, figsize=(6, 8))
    # fig11.suptitle(f'{item} Distribution of Data Mean Normalization')
    # axs = axs.flatten()
    # plt.subplots_adjust(hspace=0.5)
    # for i, lst in enumerate(dataMean):
    #     axs[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axs[i].set_title(indexes[i])
    #     axs[i].set_xlim(-2, 2)
    #     # axs[i].set_ylim(0, 40)
    # axs[-1].set_xlabel('Value')
    # fig11.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')

    bp1 = ax1.boxplot(dataMean, positions=indexes, patch_artist=True, medianprops=medianprops)
    ax1.set_xticklabels(x_values_set_ticks, fontsize=8, weight='bold')
    for i, box in enumerate(bp1['boxes']):
        if indexes[i] == 40:
            box.set(facecolor='blue')
        elif indexes[i] > 40 and indexes[i] < 90:
            box.set(facecolor='green')
        elif indexes[i] >= 90 and indexes[i] < 114:
            box.set(facecolor='orange')
        elif indexes[i] == 160:
            box.set(facecolor='purple')
        else:
            box.set(facecolor='red')
    x = np.array(indexes)
    y = np.array([np.mean(values) for values in dataMean])
    z = np.polyfit(x, y, 2)
    p = np.poly1d(z)
    # ax1.plot(x, p(x), color='magenta')
    y_pred = p(x)
    r2 = r2_score(y, y_pred)
    # Display the R-squared value on the plot
    # ax1.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax1.transAxes, fontsize=12)

    # bp2 = ax2.boxplot(dataMedian, positions=indexes, patch_artist=True, medianprops=medianprops)
    # ax2.set_xticklabels(x_values_set_ticks, fontsize=8, weight='bold')
    # for i, box in enumerate(bp2['boxes']):
    #     if indexes[i] == 40:
    #         box.set(facecolor='blue')
    #     elif indexes[i] > 40 and indexes[i] < 90:
    #         box.set(facecolor='green')
    #     elif indexes[i] >= 90 and indexes[i] < 114:
    #         box.set(facecolor='orange')
    #     elif indexes[i] == 160:
    #         box.set(facecolor='purple')
    #     else:
    #         box.set(facecolor='red')
    # x = np.array(indexes)
    # y = np.array([np.mean(values) for values in dataMedian])
    # z = np.polyfit(x, y, 2)
    # p = np.poly1d(z)
    # ax2.plot(x, p(x), color='magenta')
    # y_pred = p(x)
    # r2 = r2_score(y, y_pred)
    # # Display the R-squared value on the plot
    # ax2.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax2.transAxes, fontsize=12)
    # bp3 = ax3.boxplot(dataStd, positions=indexes, patch_artist=True, medianprops=medianprops)
    # ax3.set_xticklabels(x_values_set_ticks, fontsize=8, weight='bold')
    # for i, box in enumerate(bp3['boxes']):
    #     if indexes[i] == 40:
    #         box.set(facecolor='blue')
    #     elif indexes[i] > 40 and indexes[i] < 90:
    #         box.set(facecolor='green')
    #     elif indexes[i] >= 90 and indexes[i] < 114:
    #         box.set(facecolor='orange')
    #     elif indexes[i] == 160:
    #         box.set(facecolor='purple')
    #     else:
    #         box.set(facecolor='red')
    # x = np.array(indexes)
    # y = np.array([np.mean(values) for values in dataStd])
    # z = np.polyfit(x, y, 2)
    # p = np.poly1d(z)
    # ax3.plot(x, p(x), color='magenta')
    # y_pred = p(x)
    # r2 = r2_score(y, y_pred)
    # # Display the R-squared value on the plot
    # ax3.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax3.transAxes, fontsize=12)
    # fig.subplots_adjust(hspace=0.5)

    box_plot_data = []
    box_plot_dataNone = []
    box_plot_dataCached = []
    box_plot_dataLog = []

    indexes = [1,2,3,4,5]
    order = ['HGPS', 'HC', 'Premanifest', 'Mild', 'Severe']
    for key, value in all[item].items():
        box_plot_data.append(value)
    for key, value in allNone[item].items():
        box_plot_dataNone.append(value)
    for key, value in allCached[item].items():
        box_plot_dataCached.append(value)
    for key, value in allLog[item].items():
        box_plot_dataLog.append(value)

    fig22C, axs2C = plt.subplots(len(box_plot_dataCached), 1, figsize=(6, 8))
    plt.subplots_adjust(hspace=0.5)
    fig22C.suptitle(f'{item} Distribution of Data Raw')
    for i, lst in enumerate(box_plot_dataCached):
        axs2C[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
        axs2C[i].set_title(order[i])
        axs2C[i].set_xlim(-2, 2)
        # axs2[i].set_ylim(0, 40)
    axs2C[-1].set_xlabel('Value')
    fig22C.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    fig22C.savefig(savePlots + fr'\{item}_Raw.png', bbox_inches='tight')

    fig22N, axs2N = plt.subplots(len(box_plot_dataNone), 1, figsize=(6, 8))
    plt.subplots_adjust(hspace=0.5)
    fig22N.suptitle(f'{item} Distribution of Data Outlier Removed')
    for i, lst in enumerate(box_plot_dataNone):
        axs2N[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
        axs2N[i].set_title(order[i])
        axs2N[i].set_xlim(-2, 2)
        # axs2[i].set_ylim(0, 40)
    axs2N[-1].set_xlabel('Value')
    fig22N.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    fig22N.savefig(savePlots + fr'\{item}_OL.png', bbox_inches='tight')

    # fig22L, axs2L = plt.subplots(len(box_plot_dataLog), 1, figsize=(6, 8))
    # plt.subplots_adjust(hspace=0.5)
    # fig22L.suptitle(f'{item} Distribution of Data Log of data')
    # for i, lst in enumerate(box_plot_dataLog):
    #     axs2L[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    #     axs2L[i].set_title(order[i])
    #     axs2L[i].set_xlim(-2, 2)
    #     # axs2[i].set_ylim(0, 40)
    # axs2L[-1].set_xlabel('Value')
    # fig22L.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    # fig22L.savefig(savePlots + fr'\{item}_Log.png', bbox_inches='tight')

    fig22, axs2 = plt.subplots(len(box_plot_data), 1, figsize=(6, 8))
    plt.subplots_adjust(hspace=0.5)
    fig22.suptitle(f'{item} Distribution of Data Mean Normalization')
    for i, lst in enumerate(box_plot_data):
        axs2[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
        axs2[i].set_title(order[i])
        axs2[i].set_xlim(-2, 2)
        # axs2[i].set_ylim(0, 40)
    axs2[-1].set_xlabel('Value')
    fig22.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    fig22.savefig(savePlots + fr'\{item}_Mean.png', bbox_inches='tight')

    bp2 = ax2.boxplot(box_plot_data, positions=[5,1,2,3,4], patch_artist=True, medianprops=medianprops)
    ax2.set_xticklabels(all[item].keys())
    ax2.set_ylabel('Feature Value')
    for i, box in enumerate(bp2['boxes']):
        if indexes[i] == 2:
            box.set(facecolor='blue', alpha=0.5)
        elif indexes[i] == 3:
            box.set(facecolor='green', alpha=0.5)
        elif indexes[i] == 4:
            box.set(facecolor='orange', alpha=0.5)
        elif indexes[i] == 5:
            box.set(facecolor='red', alpha=0.5)
        else:
            box.set(facecolor='purple', alpha=0.5)
    x = np.array(indexes)
    # l = [np.mean(box_plot_data[4])], [np.mean(box_plot_data[0])], [np.mean(box_plot_data[1])], [np.mean(box_plot_data[2])], [np.mean(box_plot_data[3])]
    # y = np.array(l)
    y = np.array([np.mean(values) for values in box_plot_data])
    y = np.roll(y, -1)
    z = np.polyfit(x, y, 2)
    p = np.poly1d(z)
    # ax2.plot(x, p(x), color='magenta')
    y_pred = p(x)
    r2 = r2_score(y, y_pred)

    # Display the R-squared value on the plot
    # ax2.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax2.transAxes, fontsize=12)
    fig.savefig(savePlots + fr'\{item}_Box.png', bbox_inches='tight',  dpi=900)

    # plt.show()

    fig2, ax3 = plt.subplots(1, 1)

    bp2 = ax3.boxplot(box_plot_data, positions=[5, 1, 2, 3, 4], patch_artist=True, medianprops=medianprops)
    ax3.set_xticklabels(all[item].keys())
    current_max_y = 1.16#np.max([item.get_ydata()[1] for item in bp2["whiskers"]])  # Get the y-coordinate of the highest whisker
    offset = 0
    ax3.set_ylabel('Normalized Value')
    for i, box in enumerate(bp2['boxes']):
        if indexes[i] == 2:
            box.set(facecolor='blue', alpha=0.5)
        elif indexes[i] == 3:
            box.set(facecolor='green', alpha=0.5)
        elif indexes[i] == 4:
            box.set(facecolor='orange', alpha=0.5)
        elif indexes[i] == 5:
            box.set(facecolor='red', alpha=0.5)
        else:
            box.set(facecolor='purple', alpha=0.5)
    x = np.array(indexes)
    # l = [np.mean(box_plot_data[4])], [np.mean(box_plot_data[0])], [np.mean(box_plot_data[1])], [np.mean(box_plot_data[2])], [np.mean(box_plot_data[3])]
    # y = np.array(l)
    y = np.array([np.median(values) for values in box_plot_data])
    y = np.roll(y, -1)
    z = np.polyfit(x, y, 2)
    p = np.poly1d(z)
    # ax3.plot(x, p(x), color='magenta')
    y_pred = p(x)
    r2 = r2_score(y, y_pred)

    # Display the R-squared value on the plot
    # ax3.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax3.transAxes, fontsize=12)
    ax3.set_title('Form Factor')

    p_value_dict = extract_pValue_return(all[item])
    # Annotate p-values with stars
    for i, key1 in enumerate(all[item].keys()):
        for j, key2 in enumerate(all[item].keys()):
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
                    incr = 0.025
                    offset += current_max_y * incr  # Increment by an additional 5% of the max y-value

    fig2.savefig(savePlots + fr'\{item}_Box_Grouped.png', bbox_inches='tight', dpi=900)    # ax.legend()
    # ax.set_xlabel('Label')
    # ax.set_ylabel('Value')
    # ax.set_title(f'Scatter Plot of Stats for Item {item}')
    # plt.show()


