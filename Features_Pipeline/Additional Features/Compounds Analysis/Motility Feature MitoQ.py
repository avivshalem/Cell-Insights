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

df_migration = pd.read_excel(r'\Classifier\MigrationMitoQData.xlsx')
Dp_feature_all = df_migration.Dp
Dtot_feature_all = df_migration.Dtot
Pp_feature_all = df_migration.Pp
Psi_feature_all = df_migration.Psi
Pnp_feature_all = df_migration.Pnp
Dnp_feature_all = df_migration.Dnp
MSD_feature_all = df_migration.MSD
Sp_feature_all = df_migration.Sp
Snp_feature_all = df_migration.Snp

Dp_CAP_all = np.round(df_migration.CAP)
for i in range(len(Dp_CAP_all)):
    if Dp_CAP_all[i] < 20 and Dp_CAP_all[i] != 0:
        Dp_CAP_all[i] = np.round(df_migration.CAP[i], 1)


bins = [-1, 0, 20, 1000]  # Bins are inclusive of right value by default
labels = ['HC', 'MitoQ 1µM', 'HD']

groups = pd.cut(Dp_CAP_all, bins=bins, labels=labels)
df = pd.DataFrame({'Feature': Dtot_feature_all, 'Group': groups})
mean_group_HC = df[df['Group'] == 'HC']['Feature'].mean()
df['Normalized Feature'] = df['Feature'] / mean_group_HC

labels_ordered = ['HC', 'HD', 'MitoQ 1µM']
colors_ordered = ['blue', 'red', 'orange']

data_to_plot_ordered = [df.loc[df['Group'] == label, 'Normalized Feature'] for label in labels_ordered]




fig2, ax3 = plt.subplots(1, 1)
medianprops = dict(linestyle='-', linewidth=1.5, color='black')

bp2 = ax3.boxplot(data_to_plot_ordered, patch_artist=True, medianprops=medianprops, whis=5)
# ax3.set_xticklabels('HC', 'Premainfest', 'Mild', 'Severe', 'HGPS')
ax3.set_ylabel('Normalized Value')
for box, color in zip(bp2['boxes'], colors_ordered):
    box.set(facecolor=color, alpha=0.5)

x = np.array([1,2,3])
median_group_values = [np.median(df.loc[df['Group'] == label, 'Normalized Feature']) for label in labels_ordered]
z = np.polyfit(x, median_group_values, 2)
p = np.poly1d(z)
ax3.plot(x, p(x), color='magenta')
y_pred = p(x)
r2 = r2_score(median_group_values, y_pred)

# ax3.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax3.transAxes, fontsize=12)
ax3.set_title('Total Diffusivity')
ax3.set_xticklabels(labels_ordered)
ax3.set_ylim([-20, 10])
savePlots = r'\Motility_MitoQ'
os.makedirs(savePlots, exist_ok=True)
fig2.savefig(savePlots + fr'\Dtot_MitoQ_Box_Grouped.png', bbox_inches='tight', dpi=900)    # ax.legend()

# fig2.savefig(savePlots + fr'\Dtot_Box_Grouped.png', bbox_inches='tight', dpi=900)    # ax.legend()
# Dp_CAP = np.unique(Dp_CAP_all)
# Dp_CAP[-1] = 201.
#
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
#
# ax1.set_title(f"Box Plot Mean for Dtot")
# ax1.set_xlabel("CAP")
# ax1.set_ylabel("Feature Value")
#
# indexes = np.asarray(Dp_CAP)
# indexes[0] = 40
# indexes[-1] = 160
#
# x_values_set_ticks = [str(x) for x in Dp_CAP]
# x_values_set_ticks[0] = 'HC'
# x_values_set_ticks[-1] = 'HGPS'
#
# medianprops = dict(linestyle='-', linewidth=1.5, color='black')
#
# # bp1 = ax1.boxplot(Dtot_feature_all, positions=indexes, patch_artist=True, medianprops=medianprops)
# # ax1.set_xticklabels(x_values_set_ticks, fontsize=8, weight='bold')
# # for i, box in enumerate(bp1['boxes']):
# #     if indexes[i] == 40:
# #         box.set(facecolor='blue')
# #     elif indexes[i] > 40 and indexes[i] < 90:
# #         box.set(facecolor='green')
# #     elif indexes[i] >= 90 and indexes[i] < 114:
# #         box.set(facecolor='orange')
# #     elif indexes[i] == 160:
# #         box.set(facecolor='purple')
# #     else:
# #         box.set(facecolor='red')
# # x = np.array(indexes)
# # y = np.array([np.mean(values) for values in Dtot_feature_all])
# # z = np.polyfit(x, y, 2)
# # p = np.poly1d(z)
# # ax1.plot(x, p(x), color='magenta')
# # y_pred = p(x)
# # r2 = r2_score(y, y_pred)
# # # Display the R-squared value on the plot
# # ax1.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax1.transAxes, fontsize=12)
#
# # bp2 = ax2.boxplot(dataMedian, positions=indexes, patch_artist=True, medianprops=medianprops)
# # ax2.set_xticklabels(x_values_set_ticks, fontsize=8, weight='bold')
# # for i, box in enumerate(bp2['boxes']):
# #     if indexes[i] == 40:
# #         box.set(facecolor='blue')
# #     elif indexes[i] > 40 and indexes[i] < 90:
# #         box.set(facecolor='green')
# #     elif indexes[i] >= 90 and indexes[i] < 114:
# #         box.set(facecolor='orange')
# #     elif indexes[i] == 160:
# #         box.set(facecolor='purple')
# #     else:
# #         box.set(facecolor='red')
# # x = np.array(indexes)
# # y = np.array([np.mean(values) for values in dataMedian])
# # z = np.polyfit(x, y, 2)
# # p = np.poly1d(z)
# # ax2.plot(x, p(x), color='magenta')
# # y_pred = p(x)
# # r2 = r2_score(y, y_pred)
# # # Display the R-squared value on the plot
# # ax2.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax2.transAxes, fontsize=12)
# # bp3 = ax3.boxplot(dataStd, positions=indexes, patch_artist=True, medianprops=medianprops)
# # ax3.set_xticklabels(x_values_set_ticks, fontsize=8, weight='bold')
# # for i, box in enumerate(bp3['boxes']):
# #     if indexes[i] == 40:
# #         box.set(facecolor='blue')
# #     elif indexes[i] > 40 and indexes[i] < 90:
# #         box.set(facecolor='green')
# #     elif indexes[i] >= 90 and indexes[i] < 114:
# #         box.set(facecolor='orange')
# #     elif indexes[i] == 160:
# #         box.set(facecolor='purple')
# #     else:
# #         box.set(facecolor='red')
# # x = np.array(indexes)
# # y = np.array([np.mean(values) for values in dataStd])
# # z = np.polyfit(x, y, 2)
# # p = np.poly1d(z)
# # ax3.plot(x, p(x), color='magenta')
# # y_pred = p(x)
# # r2 = r2_score(y, y_pred)
# # # Display the R-squared value on the plot
# # ax3.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax3.transAxes, fontsize=12)
# # fig.subplots_adjust(hspace=0.5)
#
# box_plot_data = []
# box_plot_dataNone = []
# box_plot_dataCached = []
# box_plot_dataLog = []
#
# indexes = [1,2,3,4,5]
# order = ['HGPS', 'HC', 'Premanifest', 'Mild', 'Severe']
# # for key, value in all[item].items():
# #     box_plot_data.append(value)
# # for key, value in allNone[item].items():
# #     box_plot_dataNone.append(value)
# # for key, value in allCached[item].items():
# #     box_plot_dataCached.append(value)
# # for key, value in allLog[item].items():
# #     box_plot_dataLog.append(value)
#
# # fig22C, axs2C = plt.subplots(len(box_plot_dataCached), 1, figsize=(6, 8))
# # plt.subplots_adjust(hspace=0.5)
# # fig22C.suptitle(f'{item} Distribution of Data Raw')
# # for i, lst in enumerate(box_plot_dataCached):
# #     axs2C[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
# #     axs2C[i].set_title(order[i])
# #     axs2C[i].set_xlim(-2, 2)
# #     # axs2[i].set_ylim(0, 40)
# # axs2C[-1].set_xlabel('Value')
# # fig22C.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
# # fig22C.savefig(savePlots + fr'\{item}_Raw.png', bbox_inches='tight')
#
# # fig22N, axs2N = plt.subplots(len(box_plot_dataNone), 1, figsize=(6, 8))
# # plt.subplots_adjust(hspace=0.5)
# # fig22N.suptitle(f'{item} Distribution of Data Outlier Removed')
# # for i, lst in enumerate(box_plot_dataNone):
# #     axs2N[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
# #     axs2N[i].set_title(order[i])
# #     axs2N[i].set_xlim(-2, 2)
# #     # axs2[i].set_ylim(0, 40)
# # axs2N[-1].set_xlabel('Value')
# # fig22N.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
# # fig22N.savefig(savePlots + fr'\{item}_OL.png', bbox_inches='tight')
#
# # fig22L, axs2L = plt.subplots(len(box_plot_dataLog), 1, figsize=(6, 8))
# # plt.subplots_adjust(hspace=0.5)
# # fig22L.suptitle(f'{item} Distribution of Data Log of data')
# # for i, lst in enumerate(box_plot_dataLog):
# #     axs2L[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
# #     axs2L[i].set_title(order[i])
# #     axs2L[i].set_xlim(-2, 2)
# #     # axs2[i].set_ylim(0, 40)
# # axs2L[-1].set_xlabel('Value')
# # fig22L.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
# # fig22L.savefig(savePlots + fr'\{item}_Log.png', bbox_inches='tight')
#
# # fig22, axs2 = plt.subplots(len(box_plot_data), 1, figsize=(6, 8))
# # plt.subplots_adjust(hspace=0.5)
# # fig22.suptitle(f'{item} Distribution of Data Mean Normalization')
# # for i, lst in enumerate(box_plot_data):
# #     axs2[i].hist(lst, bins=len(lst), edgecolor='black', color='skyblue', alpha=1)
# #     axs2[i].set_title(order[i])
# #     axs2[i].set_xlim(-2, 2)
# #     # axs2[i].set_ylim(0, 40)
# # axs2[-1].set_xlabel('Value')
# # fig22.text(0.04, 0.5, 'Frequency', ha='center', va='center', rotation='vertical')
# # fig22.savefig(savePlots + fr'\{item}_Mean.png', bbox_inches='tight')
#
# # bp2 = ax2.boxplot(box_plot_data, positions=[5,1,2,3,4], patch_artist=True, medianprops=medianprops)
# # # ax2.set_xticklabels(all[item].keys())
# # ax2.set_ylabel('Feature Value')
# # for i, box in enumerate(bp2['boxes']):
# #     if indexes[i] == 2:
# #         box.set(facecolor='blue', alpha=0.5)
# #     elif indexes[i] == 3:
# #         box.set(facecolor='green', alpha=0.5)
# #     elif indexes[i] == 4:
# #         box.set(facecolor='orange', alpha=0.5)
# #     elif indexes[i] == 5:
# #         box.set(facecolor='red', alpha=0.5)
# #     else:
# #         box.set(facecolor='purple', alpha=0.5)
# # x = np.array(indexes)
# # # l = [np.mean(box_plot_data[4])], [np.mean(box_plot_data[0])], [np.mean(box_plot_data[1])], [np.mean(box_plot_data[2])], [np.mean(box_plot_data[3])]
# # # y = np.array(l)
# # y = np.array([np.mean(values) for values in box_plot_data])
# # y = np.roll(y, -1)
# # z = np.polyfit(x, y, 2)
# # p = np.poly1d(z)
# # ax2.plot(x, p(x), color='magenta')
# # y_pred = p(x)
# # r2 = r2_score(y, y_pred)
# #
# # # Display the R-squared value on the plot
# # ax2.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax2.transAxes, fontsize=12)
# # fig.savefig(savePlots + fr'\{item}_Box.png', bbox_inches='tight',  dpi=900)
#
# # plt.show()
#
# fig2, ax3 = plt.subplots(1, 1)
#
# bp2 = ax3.boxplot(box_plot_data, positions=[5, 1, 2, 3, 4], patch_artist=True, medianprops=medianprops)
# ax3.set_xticklabels('HC', 'Premainfest', 'Mild', 'Severe', 'HGPS')
# ax3.set_ylabel('Normalized Value')
# for i, box in enumerate(bp2['boxes']):
#     if indexes[i] == 2:
#         box.set(facecolor='blue', alpha=0.5)
#     elif indexes[i] == 3:
#         box.set(facecolor='green', alpha=0.5)
#     elif indexes[i] == 4:
#         box.set(facecolor='orange', alpha=0.5)
#     elif indexes[i] == 5:
#         box.set(facecolor='red', alpha=0.5)
#     else:
#         box.set(facecolor='purple', alpha=0.5)
# x = np.array(indexes)
# # l = [np.mean(box_plot_data[4])], [np.mean(box_plot_data[0])], [np.mean(box_plot_data[1])], [np.mean(box_plot_data[2])], [np.mean(box_plot_data[3])]
# # y = np.array(l)
# y = np.array([np.mean(values) for values in box_plot_data])
# y = np.roll(y, -1)
# z = np.polyfit(x, y, 2)
# p = np.poly1d(z)
# ax3.plot(x, p(x), color='magenta')
# y_pred = p(x)
# r2 = r2_score(y, y_pred)
#
# # Display the R-squared value on the plot
# ax3.text(0.05, 0.9, f'R-squared = {r2:.2f}', transform=ax3.transAxes, fontsize=12)
# ax3.set_title('Form Factor')
# savePlots = r'C:\Cell_Migration_MATLAB\IllustratorPlots_4\Motility'
# os.makedirs(savePlots, exist_ok=True)
# fig2.savefig(savePlots + fr'\Dtot_Box_Grouped.png', bbox_inches='tight', dpi=900)    # ax.legend()
# ax.set_xlabel('Label')
# ax.set_ylabel('Value')
# ax.set_title(f'Scatter Plot of Stats for Item {item}')
# plt.show()


