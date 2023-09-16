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

CounterXL = 0
df1 = pd.read_excel(r'\output_latB.xlsx')
df2 = pd.read_excel(r'\output_MitoQ_Severe.xlsx')
df3 = pd.read_excel(r'\output_MitoQ_Risk.xlsx')
df4 = pd.read_excel(r'\output_MitoQ.xlsx')

df = pd.concat([df1, df2, df3, df4])
df_new = pd.DataFrame()
savePlots = r'\Form_Factor_LatB'
os.makedirs(savePlots, exist_ok=True)
cols = df.iloc[:, 0]
new_cols = []
for itemIndex, item in enumerate(cols):
    if 'HCMedia' in item or 'HDMedia' in item or 'MitoQ' in item:
        continue
    new_cols.append(item)
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
    if 'HCMedia' in new_cols[col]:
        continue
        type_cols.append('HCMedia')
    elif 'HC' in new_cols[col]:
        type_cols.append('HC')
    elif 'HDMedia' in new_cols[col]:
        continue
        type_cols.append('HDMedia')
    elif 'HD' in new_cols[col]:
        type_cols.append('HD')
    elif 'LatB' in new_cols[col]:
        type_cols.append('LatB 60nM')
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
        # outliers = model.predict(data.reshape(-1, 1))    bp2 = ax2.boxplot(box_plot_data, positions=[1,2,3,4,5], patch_artist=True, medianprops=medianprops)
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

    fig, ax2 = plt.subplots(1, 1)
    # plt.subplots_adjust(hspace=0.5)
    #
    # ax1.set_title(f"Box Plot Mean for {item}")
    # ax1.set_xlabel("CAP")
    # ax1.set_ylabel("Feature Value")
    # # fig1, (ax11, ax22) = plt.subplots(2, 1, figsize=(6, 8))
    # # ax11.set_title(f"Distribution plot for {item} by CAP")
    # # ax11.set_xlabel("CAP")
    # # ax11.set_ylabel("Occurrences")
    # # ax22.set_title(f"Distribution plot for {item}")
    # # ax22.set_ylabel("Occurrences")
    #
    # # x_valuesMedian = list(dicMedian.keys())
    # # ax2.set_title(f"Box Plot Median for {item}")
    # # ax2.set_xlabel("X values")
    # # ax2.set_ylabel("Y values")
    # # x_valuesStd = list(dicStd.keys())
    # # ax3.set_title(f"Box Plot STD for {item}")
    # # ax3.set_xlabel("X values")
    # # ax3.set_ylabel("Y values")
    #
    # # Loop through the x_values and plot the corresponding y_values
    # # for x_value in x_valuesMean:
    # #     ax1.scatter([x_value] * len(dicMean[x_value]), dicMean[x_value])
    # # for x_value in x_valuesMedian:
    # #     ax2.scatter([x_value] * len(dicMedian[x_value]), dicMedian[x_value])
    # # for x_value in x_valuesStd:
    # #     ax3.scatter([x_value] * len(dicStd[x_value]), dicStd[x_value])
    # # fig.subplots_adjust(hspace=0.5)
    # # plt.show()
    # x_values_set = natsort.natsorted(list(dicMean.keys()))
    # # x_values_set = ['HC', 'Premanifest', 'Mild', 'Severe']
    # float_list = [float(s) for s in x_values_set]
    # indexes = np.asarray(float_list)
    # indexes[0] = 40
    # indexes[-1] = 160
    #
    # x_values_set_ticks = x_values_set.copy()
    # x_values_set_ticks[0] = 'HC'
    # x_values_set_ticks[-1] = 'HGPS'
    #
    # dataMean, dataMeanNone, dataMeanCached, dataMeanLog, dataMedian, dataStd = [], [], [], [], [], []
    # for x_value in x_values_set:
    #     currValue = dicMean[x_value]
    #     dataMeanNone.append(currValue)
    #     # currValue = np.log(currValue)
    #     dataMeanLog.append(np.log(dicMean[x_value]))
    #     # res = [((x - mean[item]) / std[item]) for x in currValue]
    #     res = [(x / mean[item]) for x in currValue]
    #     dataMeanCached.append(dicMeanCached[x_value])
    #     dataMean.append(res)
    #     # ax1.scatter([x_value] * len(dicMean[x_value]), dicMean[x_value])
    #     # ax1.scatter(x_value, dicMean[x_value])
    # # for x_value in x_values_set:
    # #     dataMedian.append(dicMedian[x_value])
    # #     # ax2.scatter([x_value] * len(dicMedian[x_value]), dicMedian[x_value])
    # #     # ax2.scatter(x_value, dicMedian[x_value])
    # # for x_value in x_values_set:
    # #     dataStd.append(dicStd[x_value])
    # #     # ax3.scatter([x_value] * len(dicStd[x_value]), dicStd[x_value])
    # #     # ax3.scatter(x_value, dicStd[x_value])
    # # for i, group in enumerate(dataMean):
    # #     # Combine all the values into a single list
    # #     values = [val for val in group]
    # #     sns.displot(values, ax=ax11)
    #
    # # fig11C, axsC = plt.subplots(int(len(dataMeanCached) / 4), 4, figsize=(6, 8))
    # # plt.subplots_adjust(hspace=0.5)
    # # fig11C.suptitle(f'{item} Distribution of Data Raw')
    # # axsC = axsC.flatten()
    # # for i, lst in enumerate(dataMeanCached):
    # #     axsC[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    # #     axsC[i].set_title(indexes[i])
    # #     axsC[i].set_xlim(-2, 2)
    # #     # axs2[i].set_ylim(0, 40)
    # # axsC[-1].set_xlabel('Value')
    # # fig11C.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    # #
    # # fig11N, axsN = plt.subplots(int(len(dataMeanNone)/4), 4, figsize=(6, 8))
    # # plt.subplots_adjust(hspace=0.5)
    # # fig11N.suptitle(f'{item} Distribution of Data Outlier Removed')
    # # axsN = axsN.flatten()
    # # for i, lst in enumerate(dataMeanNone):
    # #     axsN[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    # #     axsN[i].set_title(indexes[i])
    # #     axsN[i].set_xlim(-2, 2)
    # #     # axs2[i].set_ylim(0, 40)
    # # axsN[-1].set_xlabel('Value')
    # # fig11N.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    # #
    # # fig11L, axsL = plt.subplots(int(len(dataMeanLog) / 4), 4, figsize=(6, 8))
    # # plt.subplots_adjust(hspace=0.5)
    # # fig11L.suptitle(f'{item} Distribution of Log of Data')
    # # axsL = axsL.flatten()
    # # for i, lst in enumerate(dataMeanLog):
    # #     axsL[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    # #     axsL[i].set_title(indexes[i])
    # #     axsL[i].set_xlim(-2, 2)
    # #     # axs2[i].set_ylim(0, 40)
    # # axsL[-1].set_xlabel('Value')
    # # fig11L.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    # #
    # # fig11, axs = plt.subplots(int(len(dataMean)/4), 4, figsize=(6, 8))
    # # fig11.suptitle(f'{item} Distribution of Data Mean Normalization')
    # # axs = axs.flatten()
    # # plt.subplots_adjust(hspace=0.5)
    # # for i, lst in enumerate(dataMean):
    # #     axs[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
    # #     axs[i].set_title(indexes[i])
    # #     axs[i].set_xlim(-2, 2)
    # #     # axs[i].set_ylim(0, 40)
    # # axs[-1].set_xlabel('Value')
    # # fig11.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
    #
    # bp1 = ax1.boxplot(dataMean, positions=indexes, patch_artist=True, medianprops=medianprops)
    # ax1.set_xticklabels(x_values_set_ticks, fontsize=8, weight='bold')
    # for i, box in enumerate(bp1['boxes']):
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
    # y = np.array([np.mean(values) for values in dataMean])
    # z = np.polyfit(x, y, 2)
    # p = np.poly1d(z)
    # ax1.plot(x, p(x), color='magenta')
    # y_pred = p(x)
    # r2 = r2_score(y, y_pred)
    # # Display the R-squared value on the plot
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

    indexes = [1,2,3]
    order = ['HC', 'HD','LatB 60nM']
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

    bp2 = ax2.boxplot(box_plot_data, positions=[1,2,3], patch_artist=True, medianprops=medianprops)
    ax2.set_xticklabels(all[item].keys())
    ax2.set_ylabel('Feature Value')
    for i, box in enumerate(bp2['boxes']):
        if indexes[i] == 1:
            box.set(facecolor='blue', alpha=0.5)
        elif indexes[i] == 2:
            box.set(facecolor='red', alpha=0.5)
        elif indexes[i] == 3:
            box.set(facecolor='orange', alpha=0.5)
        elif indexes[i] == 5:
            box.set(facecolor='red', alpha=0.5)
        else:
            box.set(facecolor='purple', alpha=0.5)
    x = np.array(indexes)
    # l = [np.mean(box_plot_data[4])], [np.mean(box_plot_data[0])], [np.mean(box_plot_data[1])], [np.mean(box_plot_data[2])], [np.mean(box_plot_data[3])]
    # y = np.array(l)
    y = np.array([np.mean(values) for values in box_plot_data])
    # y = np.roll(y, -1)
    z = np.polyfit(x, y, 2)
    p = np.poly1d(z)
    ax2.plot(x, p(x), color='magenta')
    y_pred = p(x)
    r2 = r2_score(y, y_pred)

    # Display the R-squared value on the plot
    # ax2.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax2.transAxes, fontsize=12)
    ax2.set_title(f"Box Plot Mean for {item}")
    fig.savefig(savePlots + fr'\{item}_Box.png', bbox_inches='tight', dpi=900)

    # plt.show()

    # ax.legend()
    # ax.set_xlabel('Label')
    # ax.set_ylabel('Value')
    # ax.set_title(f'Scatter Plot of Stats for Item {item}')
    # plt.show()


