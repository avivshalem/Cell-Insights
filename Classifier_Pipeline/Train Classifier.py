import os
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from keras.models import Sequential
from keras.layers import Dense, Dropout
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
from itertools import chain, combinations
from keras.callbacks import ModelCheckpoint#

# EXAMPLE DATA ONLY:
Dp_CAP = np.load(r'example_data\Dp_CAP.npy')
Dp_feature_Train = np.load(r'example_data\Dp_feature_Train.npy')
ratio_width_to_length_mean_feature_train = np.load(r'example_data\ratio_width_to_length_mean_feature_train.npy')
total_area_sum_feature_train = np.load(r'example_data\total_area_sum_feature_train.npy')
form_factor_feature_train = np.load(r'example_data\form_factor_feature_train.npy')
Dp_feature_Test = np.load(r'example_data\Dp_feature_Test.npy')
ratio_width_to_length_mean_feature_test = np.load(r'example_data\ratio_width_to_length_mean_feature_test.npy')
total_area_sum_feature_test = np.load(r'example_data\total_area_sum_feature_test.npy')
form_factor_feature_test = np.load(r'example_data\form_factor_feature_test.npy')
Dtot_feature_Train = np.load(r'example_data\Dtot_feature_Train.npy')
Pp_feature_Train = np.load(r'example_data\Pp_feature_Train.npy')
Psi_feature_Train = np.load(r'example_data\Psi_feature_Train.npy')
Pnp_feature_Train = np.load(r'example_data\Pnp_feature_Train.npy')
Dnp_feature_Train = np.load(r'example_data\Dnp_feature_Train.npy')
MSD_feature_Train = np.load(r'example_data\MSD_feature_Train.npy')
Sp_feature_Train = np.load(r'example_data\Sp_feature_Train.npy')
Snp_feature_Train = np.load(r'example_data\Snp_feature_Train.npy')
Dtot_feature_Test = np.load(r'example_data\Dtot_feature_Test.npy')
Pp_feature_Test = np.load(r'example_data\Pp_feature_Test.npy')
Psi_feature_Test = np.load(r'example_data\Psi_feature_Test.npy')
MSD_feature_Test = np.load(r'example_data\MSD_feature_Test.npy')
Pnp_feature_Test = np.load(r'example_data\Pnp_feature_Test.npy')
Dnp_feature_Test = np.load(r'example_data\Dnp_feature_Test.npy')
Snp_feature_Test = np.load(r'example_data\Snp_feature_Test.npy')
Sp_feature_Test = np.load(r'example_data\Sp_feature_Test.npy')
# FINISH EXAMPLE DATA READ

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

runRegression = False
if runRegression:
    # Fit and evaluate Polynomial Regression
    poly = PolynomialFeatures(degree=2)
    poly_train_data = poly.fit_transform(train_data_scaled)
    poly_test_data = poly.transform(test_data_scaled)  # apply same transformation to test data
    model_poly = LinearRegression()
    model_poly.fit(poly_train_data, Dp_CAP)
    predictions_poly = model_poly.predict(poly_test_data)  # make predictions on the test data

    # Fit and evaluate Decision Tree Regression
    model_tree = DecisionTreeRegressor()
    model_tree.fit(train_data_scaled, Dp_CAP)
    predictions_tree = model_tree.predict(test_data_scaled)  # make predictions on the test data

    # Fit and evaluate SVR
    model_svr = SVR()
    model_svr.fit(train_data_scaled, Dp_CAP)
    predictions_svr = model_svr.predict(test_data_scaled)  # make predictions on the test data

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


runClassifier = False
if runClassifier:
    # Fit and evaluate Decision Tree Classifier
    model_tree_classifier = DecisionTreeClassifier()
    model_tree_classifier.fit(train_data_scaled, Dp_CAP_Discrete)
    predictions_tree_classifier = model_tree_classifier.predict(test_data_scaled)  # make predictions on the test data

    # Fit and evaluate Random Forest Classifier
    model_random_forest = RandomForestClassifier()
    model_random_forest.fit(train_data_scaled, Dp_CAP_Discrete)
    predictions_random_forest = model_random_forest.predict(test_data_scaled)  # make predictions on the test data

    # Fit and evaluate SVC
    model_svc = SVC()
    model_svc.fit(train_data_scaled, Dp_CAP_Discrete)
    predictions_svc = model_svc.predict(test_data_scaled)  # make predictions on the test data

runNetwork = False
if runNetwork:
    path = r'\Models_Classifier'
    os.makedirs(path, exist_ok=True)

    features_train = [Dp_feature_Train, ratio_width_to_length_mean_feature_train,
                      total_area_sum_feature_train, form_factor_feature_train,
                      Dtot_feature_Train, Pp_feature_Train, Psi_feature_Train,
                      Pnp_feature_Train, Dnp_feature_Train, MSD_feature_Train,
                      Sp_feature_Train, Snp_feature_Train]

    features_test = [Dp_feature_Test, ratio_width_to_length_mean_feature_test,
                     total_area_sum_feature_test, form_factor_feature_test,
                     Dtot_feature_Test, Pp_feature_Test, Psi_feature_Test,
                     Pnp_feature_Test, Dnp_feature_Test, MSD_feature_Test,
                     Sp_feature_Test, Snp_feature_Test]

    # First, we train a model with all the features
    train_data = np.concatenate(features_train, axis=1)
    test_data = np.concatenate(features_test, axis=1)

    scaler = MinMaxScaler()
    train_data_scaled = scaler.fit_transform(train_data)
    test_data_scaled = scaler.transform(test_data)

    le = LabelEncoder()
    labels = le.fit_transform(Dp_CAP_Discrete)
    # Define the model
    model = Sequential()
    model.add(Dense(256, input_dim=train_data_scaled.shape[1], activation='relu'))
    model.add(Dropout(0.5))  # dropout layer
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.5))  # dropout layer
    model.add(Dense(64, activation='relu'))
    model.add(Dropout(0.5))  # dropout layer
    model.add(Dense(32, activation='relu'))
    model.add(Dropout(0.5))  # dropout layer
    model.add(Dense(len(np.unique(labels)),
                    activation='softmax'))  # The number of nodes in the output layer should be equal to the number of classes

    # Compile the model
    model.compile(loss='sparse_categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

    history = model.fit(train_data_scaled, labels, validation_data=(test_data_scaled, labels), epochs=1000,
                        batch_size=1)

    # Save the model
    best_accuracy = max(history.history['val_accuracy'])
    model.save(rf'{path}\model_all_features_Accuracy_{best_accuracy}.h5')

    # Plot training & validation accuracy values
    plt.figure(figsize=(12, 6))
    plt.plot(history.history['accuracy'])
    plt.plot(history.history['val_accuracy'])
    plt.title('Model accuracy')
    plt.ylabel('Accuracy')
    plt.xlabel('Epoch')
    plt.yticks(np.arange(0, 1.05, 0.05))  # this line controls the y-axis
    plt.legend(['Train', 'Test'], loc='upper left')
    plt.savefig(rf'{path}\model_all_features_accuracy.png')
    plt.close()

    # Plot training & validation loss values
    plt.figure(figsize=(12, 6))
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('Model loss')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
    plt.savefig(rf'{path}\model_all_features_loss.png')
    plt.close()

    # Now, we train a model for each combination of features without one of them
    for i in range(len(features_train)):
        train_data = np.concatenate([features_train[j] for j in range(len(features_train)) if j != i], axis=1)
        test_data = np.concatenate([features_test[j] for j in range(len(features_test)) if j != i], axis=1)

        scaler = MinMaxScaler()
        train_data_scaled = scaler.fit_transform(train_data)
        test_data_scaled = scaler.transform(test_data)

        le = LabelEncoder()
        labels = le.fit_transform(Dp_CAP_Discrete)
        # Define the model
        model = Sequential()
        model.add(Dense(256, input_dim=train_data_scaled.shape[1], activation='relu'))
        model.add(Dropout(0.5))  # dropout layer
        model.add(Dense(128, activation='relu'))
        model.add(Dropout(0.5))  # dropout layer
        model.add(Dense(64, activation='relu'))
        model.add(Dropout(0.5))  # dropout layer
        model.add(Dense(32, activation='relu'))
        model.add(Dropout(0.5))  # dropout layer
        model.add(Dense(len(np.unique(labels)),
                        activation='softmax'))  # The number of nodes in the output layer should be equal to the number of classes

        # Compile the model
        model.compile(loss='sparse_categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

        history = model.fit(train_data_scaled, labels, validation_data=(test_data_scaled, labels), epochs=1000,
                            batch_size=1)

        # Save the model
        best_accuracy = max(history.history['val_accuracy'])
        model.save(rf'{path}\model_without_feature_{i}_Accuracy_{best_accuracy}.h5')

        # Plot training & validation accuracy values
        plt.figure(figsize=(12, 6))
        plt.plot(history.history['accuracy'])
        plt.plot(history.history['val_accuracy'])
        plt.title('Model accuracy')
        plt.ylabel('Accuracy')
        plt.xlabel('Epoch')
        plt.yticks(np.arange(0, 1.05, 0.05))  # this line controls the y-axis
        plt.legend(['Train', 'Test'], loc='upper left')
        plt.savefig(rf'{path}\model_without_feature_{i}_accuracy.png')
        plt.close()

        # Plot training & validation loss values
        plt.figure(figsize=(12, 6))
        plt.plot(history.history['loss'])
        plt.plot(history.history['val_loss'])
        plt.title('Model loss')
        plt.ylabel('Loss')
        plt.xlabel('Epoch')
        plt.legend(['Train', 'Test'], loc='upper left')
        plt.savefig(rf'{path}\model_without_feature_{i}_loss.png')
        plt.close()

runComb = True
if runComb:
    path = r'\Models_Classifier'
    os.makedirs(path, exist_ok=True)

    features_train = [Dp_feature_Train, ratio_width_to_length_mean_feature_train,
                      total_area_sum_feature_train, form_factor_feature_train,
                      Dtot_feature_Train, Dnp_feature_Train, MSD_feature_Train]

    features_test = [Dp_feature_Test, ratio_width_to_length_mean_feature_test,
                     total_area_sum_feature_test, form_factor_feature_test,
                     Dtot_feature_Test, Dnp_feature_Test, MSD_feature_Test]

    # Get all combinations of features
    combs = chain(*map(lambda x: combinations(range(len(features_train)), x), range(1, len(features_train) + 1)))

    for comb in combs:
        train_data = np.concatenate([features_train[i] for i in comb], axis=1)
        test_data = np.concatenate([features_test[i] for i in comb], axis=1)

        scaler = MinMaxScaler()
        train_data_scaled = scaler.fit_transform(train_data)
        test_data_scaled = scaler.transform(test_data)

        le = LabelEncoder()
        labels = le.fit_transform(Dp_CAP_Discrete)
        # Define the model
        model = Sequential()
        model.add(Dense(256, input_dim=train_data_scaled.shape[1], activation='relu'))
        model.add(Dropout(0.5))  # dropout layer
        model.add(Dense(128, activation='relu'))
        model.add(Dropout(0.5))  # dropout layer
        model.add(Dense(64, activation='relu'))
        model.add(Dropout(0.5))  # dropout layer
        model.add(Dense(32, activation='relu'))
        model.add(Dropout(0.5))  # dropout layer
        model.add(Dense(len(np.unique(labels)),
                        activation='softmax'))  # The number of nodes in the output layer should be equal to the number of classes

        # Compile the model
        model.compile(loss='sparse_categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

        # create a checkpoint object
        checkpoint = ModelCheckpoint(rf'{path}\model_{comb}.h5', monitor='val_accuracy', verbose=1,
                                     save_best_only=True, mode='max')
        callbacks_list = [checkpoint]
        history = model.fit(train_data_scaled, labels, validation_data=(test_data_scaled, labels), epochs=1000,
                            batch_size=1, callbacks=callbacks_list)

        # Save the model
        best_accuracy = max(history.history['val_accuracy'])
        # model.save(rf'{path}\model_{comb}_Accuracy_{best_accuracy}.h5')

        # Plot training & validation accuracy values
        plt.figure(figsize=(12, 6))
        plt.plot(history.history['accuracy'])
        plt.plot(history.history['val_accuracy'])
        plt.title('Model accuracy')
        plt.ylabel('Accuracy')
        plt.xlabel('Epoch')
        plt.yticks(np.arange(0, 1.05, 0.05))  # this line controls the y-axis
        plt.legend(['Train', 'Test'], loc='upper left')
        plt.savefig(rf'{path}\model_{comb}_accuracy_{best_accuracy}.png')
        plt.close()

        # Plot training & validation loss values
        plt.figure(figsize=(12, 6))
        plt.plot(history.history['loss'])
        plt.plot(history.history['val_loss'])
        plt.title('Model loss')
        plt.ylabel('Loss')
        plt.xlabel('Epoch')
        plt.legend(['Train', 'Test'], loc='upper left')
        plt.savefig(rf'{path}\model_{comb}_loss.png')
        plt.close()
