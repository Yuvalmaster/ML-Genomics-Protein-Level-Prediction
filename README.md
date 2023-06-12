# Protein Level Predictor for Synthetic Operons 🪱 🔬 🧪
This project aims to train a model to determine protein levels in a library of variants of Synthetic operons (consists of 12,900 variants). The algorithm predicts the protein levels for each variant by dividing them into bins. The predictor utilizes a set of features extracted from the operon sequences to make accurate predictions.

## Project Execution Instructions
To successfully execute the code for the protein level predictor, please follow the instructions below:

* Open the file named 'Correlation_and_predictions.m' in MATLAB.
* Inside the folder, you will find the extracted features and a table containing the features for both the known and unknown data. By default, the variable 'create' is set to 0 in the file. If you want to recreate the features tables, change the value of 'create' to 1 instead of 0.
* Press the 'RUN' button in the editor section of MATLAB to execute the code.
* After the code execution, you will have access to two important variables:
'Model_prediction_table': This variable contains the accepted protein levels for all models.
'Final_Result_prediction': This column vector contains the predicted protein levels for the unknown data. The predictions are obtained by choosing the Narrow Neural Network model. Each row in the vector corresponds to a bin, and the order of the bins is the same as specified in the file 'unknown_data_set.xlsx'.

## Loading the Results
To load the results obtained from the protein level predictor, you have two options:

* Open the file 'Final_Result_prediction.mat', which contains a column vector ('Final_Result_prediction') with the predicted protein levels for the unknown data set.
* Load the Excel file named 'Unknown_data_set_w_AverageBinPA_predictions.xlsx'. Column 28 in this file contains the predicted protein levels (AverageBinPA_Prediction), similar to the known_data_set. Please note that the Excel file contains all the features used for the prediction analysis, providing comprehensive information about the predictions.

If you have any questions or need additional assistance, please contact me for further clarification.
