## Project Description
This project aims to build a predictor that can determine protein levels in a library of variants of Synthetic operons. The library consists of 12,900 variants, and the algorithm predicts the protein levels over a sequence of each variant and divides them into bins.

## Instructions
Follow the instructions below to successfully execute the code:

Open the file named 'Correlation_and_predictions.m'.
The folder already contains the extraction of the features, and a table with the features for both the known and unknown data already exists. In this case, the variable 'create' is set to 0. If you want to recreate the features tables, change the variable 'create' in the file to 1 instead of 0.
Press 'RUN' in the editor section of MATLAB to execute the code.
The variable 'Model_prediction_table' contains the accepted protein levels for all models.
The variable 'Final_Result_prediction' is a column vector that contains the protein levels for the unknown data, obtained by choosing the Narrow Neural Network model. Each row corresponds to a bin, and the order of the bins is the same as the one specified in the file 'unknown_data_set.xlsx'.
To load the results, you have two options:
Open the file 'Final_Result_prediction.mat', which contains a column vector with the predicted protein levels for the unknown data set.
Load 'Unknown_data_set_w_AverageBinPA_predictions.xlsx'. Column 28 in this file contains the predicted protein levels (AverageBinPA_Prediction), similar to the known_data_set. Please note that the Excel file contains all the features used for the prediction analysis.

