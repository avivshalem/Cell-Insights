import os

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Define the classes and the test set inputs
classes = ["HC", "Premanifest", "Mild", "Severe", "HGPS"]
test_inputs = ["HC", "Severe MitoQ 1µM", "Severe MitoQ 1µM", "Severe MitoQ 1µM", "Severe", "Severe", "Severe"]

# Initialize a matrix with zeros
data = np.zeros((len(classes), len(test_inputs)))

# Manually set the values where the model made a prediction
# Replace these indices with the correct ones according to your data
data[0, 0] = 1  # The model predicted the first test input to be the first class
data[2, 1] = 0.75  # The model predicted the second test input to be the fourth class
data[1, 2] = 0.75  # The model predicted the first test input to be the first class
data[1, 3] = 0.75  # The model predicted the second test input to be the fourth class
data[3, 4] = 1  # The model predicted the first test input to be the first class
data[3, 5] = 1  # The model predicted the second test input to be the fourth class
data[3, 6] = 1  # The model predicted the second test input to be the fourth class

# Continue this for all predictions

# Convert the matrix to a DataFrame for easier plotting
data_df = pd.DataFrame(data, index=classes, columns=test_inputs)

# Plot the heatmap
plt.figure(figsize=(10, 7))
sns.heatmap(data_df, annot=False, cmap='Greens', cbar=False, linecolor='black', linewidths=1)
plt.xlabel("MitoQ 1µM Experiment Inputs")
plt.ylabel("Classes")
plt.xticks(rotation=30)  # This line rotates the x-axis labels by 45 degrees
plt.title('Model Results on MitoQ 1µM Experiment Data')
# plt.show()
savePlots = r'C:\Users\avivs\PycharmProjects\CellInsights\Example_Data\Cell_Insights_Large_Files\Features\Features\Classifier'
os.makedirs(savePlots, exist_ok=True)
plt.savefig(savePlots + fr'\Heatmap2.png', bbox_inches='tight', dpi=900)
