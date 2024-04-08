import os
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from keras.models import load_model
from sklearn.preprocessing import LabelEncoder
from itertools import chain, combinations
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt


# EXAMPLE DATA ONLY:
Dp_CAP = np.load(r'example_data_inference\Dp_CAP.npy')
Dp_feature_Train = np.load(r'example_data_inference\Dp_feature_Train.npy')
ratio_width_to_length_mean_feature_train = np.load(r'example_data_inference\ratio_width_to_length_mean_feature_train.npy')
total_area_sum_feature_train = np.load(r'example_data_inference\total_area_sum_feature_train.npy')
form_factor_feature_train = np.load(r'example_data_inference\form_factor_feature_train.npy')
Dp_feature_Test = np.load(r'example_data_inference\Dp_feature_Test.npy')
ratio_width_to_length_mean_feature_test = np.load(r'example_data_inference\ratio_width_to_length_mean_feature_test.npy')
total_area_sum_feature_test = np.load(r'example_data_inference\total_area_sum_feature_test.npy')
form_factor_feature_test = np.load(r'example_data_inference\form_factor_feature_test.npy')
Dtot_feature_Train = np.load(r'example_data_inference\Dtot_feature_Train.npy')
Pp_feature_Train = np.load(r'example_data_inference\Pp_feature_Train.npy')
Psi_feature_Train = np.load(r'example_data_inference\Psi_feature_Train.npy')
Pnp_feature_Train = np.load(r'example_data_inference\Pnp_feature_Train.npy')
Dnp_feature_Train = np.load(r'example_data_inference\Dnp_feature_Train.npy')
MSD_feature_Train = np.load(r'example_data_inference\MSD_feature_Train.npy')
Sp_feature_Train = np.load(r'example_data_inference\Sp_feature_Train.npy')
Snp_feature_Train = np.load(r'example_data_inference\Snp_feature_Train.npy')
Dtot_feature_Test = np.load(r'example_data_inference\Dtot_feature_Test.npy')
Pp_feature_Test = np.load(r'example_data_inference\Pp_feature_Test.npy')
Psi_feature_Test = np.load(r'example_data_inference\Psi_feature_Test.npy')
MSD_feature_Test = np.load(r'example_data_inference\MSD_feature_Test.npy')
Pnp_feature_Test = np.load(r'example_data_inference\Pnp_feature_Test.npy')
Dnp_feature_Test = np.load(r'example_data_inference\Dnp_feature_Test.npy')
Snp_feature_Test = np.load(r'example_data_inference\Snp_feature_Test.npy')
Sp_feature_Test = np.load(r'example_data_inference\Sp_feature_Test.npy')
# FINISH EXAMPLE DATA READ

fpr_boot = dict()
tpr_boot = dict()
interp_tpr_boot = dict()
roc_auc_boot = dict()
total_bootstrapped = 10000
for boot in range(total_bootstrapped):

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

    path = r'\Models_CLassifier'
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
        r'Models_Classifier\model_(1, 2, 3, 4).h5')
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


interp_tpr_boot = {k: np.mean(v, axis=0) for k, v in interp_tpr_boot.items()}
roc_auc_boot = {k: v / total_bootstrapped for k, v in roc_auc_boot.items()}

# Plot all ROC curves
plt.figure()
colors = [(66 / 255, 133 / 255, 244 / 255), (51 / 255, 0, 114 / 255),
            (244 / 255, 180 / 255, 0), (15 / 255, 157 / 255, 88 / 255), (219 / 255, 68 / 255, 55 / 255)]

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
plt.savefig(r'example_data_inference\ROC_Bootstrapped10k.png', bbox_inches='tight', dpi=900)
