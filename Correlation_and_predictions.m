clc; clear all; close all;

ModelsDir = './Models'      ; addpath(ModelsDir);
FeaturesDir = './Features'  ; addpath(FeaturesDir);
%% Find all features & Create data tables
% create feature flag
create = 0;             % Change to 1 if features tables do not exist ; change to 0 if features tables do exist

if create == 1
    [known_with_features, known_features, unknown_with_features, unknown_features] = Run_features();
else
    known_with_features     = readtable('known_with_features.xlsx','VariableNamingRule','preserve');
    known_features          = readtable('known_features.xlsx','VariableNamingRule','preserve');
    unknown_with_features   = readtable('unknown_with_features.xlsx','VariableNamingRule','preserve');
    unknown_features        = readtable('unknown_features.xlsx','VariableNamingRule','preserve'); clc;
    
end
Average_bin_PA = known_features{:,1};

%% Train Models 
% --- All models were built by Linear Regression learner App of MATLAB.
% --- The models were trained by the App, and 10-fold cross-validation.
% --- In addition for each model, and RMSE validation value is extracted to
% --- verify good validation of the model.
       
   
% SVM regression model
[SVM_Cubic_model, SVM_Cubic_validationRMSE]                                     = trainSVMCubicModel(known_features);

% Linear regression models
[Linear_Regression_Model, Linear_Regression_validationRMSE]                     = trainLinearRegressionModel(known_features);
[Robust_Linear_Regression_model, Robust_Linear_Regression_validationRMSE]       = trainRobustLinearRegressionModel(known_features);

% Neural Network models
[Narrow_Neural_Network_Model, Narrow_Neural_Network_validationRMSE]             = trainNarrowNeuralNetworkModel(known_features);
[Medium_Neural_Network_Model, Medium_Neural_Network_validationRMSE]             = trainMediumNeuralNetworkModel(known_features);

% Gaussian regression model
[Gaussian_Process_Regression_Model, Gaussian_Process_Regression_validationRMSE] = trainGaussianProcessRegressionModel(known_features);

% Ensemble Trees models
[Ensemble_Boosted_Trees_Model, Ensemble_Boosted_Trees_validationRMSE]           = trainEnsembleBoostedTreesModel(known_features);
[Ensemble_Bagged_Trees_Model, Ensemble_Bagged_Trees_validationRMSE]             = trainEnsembleBaggedTreesModel(known_features);
clc;

%% Predict models based on known data set
yfit_SVM      = SVM_Cubic_model.predictFcn(known_features(:,2:end));

yfit_Linear1  = Linear_Regression_Model.predictFcn(known_features(:,2:end));
yfit_Linear2  = Robust_Linear_Regression_model.predictFcn(known_features(:,2:end));

yfit_NN1      = Narrow_Neural_Network_Model.predictFcn(known_features(:,2:end));
yfit_NN2      = Medium_Neural_Network_Model.predictFcn(known_features(:,2:end));

yfit_Gaussian = Gaussian_Process_Regression_Model.predictFcn(known_features(:,2:end));

yfit_ET1      = Ensemble_Boosted_Trees_Model.predictFcn(known_features(:,2:end));
yfit_ET2      = Ensemble_Bagged_Trees_Model.predictFcn(known_features(:,2:end));

% Plot prediction of each model
figure; plot(Average_bin_PA,'--'); hold on;
        plot(yfit_SVM);
        plot(yfit_Linear1);
        plot(yfit_Linear2);
        plot(yfit_NN1);
        plot(yfit_NN2);
        plot(yfit_Gaussian);
        plot(yfit_ET1);
        plot(yfit_ET2);
title('Comparison between models predictions of Average Bin PA');
ylabel('Average Bin PA'); 
legend('Average Bin PA', 'Cubic SVM', 'Linear Regression', ...
       'Robust Linear Regression'   , 'Narrow Neural Network',...
       'Medium Neural Network'      , 'Gaussian Process Regression',...
       'Ensemble Boosted Trees'     , 'Ensemble Bagged Trees');

%% Spearman correlation between trained model to response
correlations = [corr(Average_bin_PA, yfit_SVM,          type = 'Spearman')
                corr(Average_bin_PA, yfit_Linear1,      type = 'Spearman')
                corr(Average_bin_PA, yfit_Linear2,      type = 'Spearman')
                corr(Average_bin_PA, yfit_NN1,          type = 'Spearman')
                corr(Average_bin_PA, yfit_NN2,          type = 'Spearman')
                corr(Average_bin_PA, yfit_Gaussian,     type = 'Spearman')
                corr(Average_bin_PA, yfit_ET1,          type = 'Spearman')
                corr(Average_bin_PA, yfit_ET2,          type = 'Spearman')];
Model_names = {'Cubic SVM'                ; 'Linear Regression'; ...
               'Robust Linear Regression' ; 'Narrow Neural Network';...
               'Medium Neural Network'    ; 'Gaussian Process Regression';...
               'Ensemble Boosted Trees'   ; 'Ensemble Bagged Trees'};

Model_prediction_correlation_table = ...
    array2table(correlations','variableNames', Model_names);

disp(' ')
disp('-------------------------------------------------------------------------------------------------------------------------------')
disp(' '); disp('<strong>Correlation Table between models </strong>')
disp(' ')
disp('-------------------------------------------------------------------------------------------------------------------------------')
disp(Model_prediction_correlation_table);


%% Unkown data predictions

Unkown_Predictions_Model_data =[
    SVM_Cubic_model.predictFcn(unknown_features)';

    Linear_Regression_Model.predictFcn(unknown_features)';
    Robust_Linear_Regression_model.predictFcn(unknown_features)';

    Narrow_Neural_Network_Model.predictFcn(unknown_features)';
    Medium_Neural_Network_Model.predictFcn(unknown_features)';

    Gaussian_Process_Regression_Model.predictFcn(unknown_features)';

    Ensemble_Boosted_Trees_Model.predictFcn(unknown_features)';
    Ensemble_Bagged_Trees_Model.predictFcn(unknown_features)']';

Model_prediction_table = ...
    array2table(Unkown_Predictions_Model_data,'variableNames', Model_names);

Final_Result_prediction = Unkown_Predictions_Model_data(:,4);
save('Final_Result_prediction.mat','Final_Result_prediction');
%% Create table

Unknown_data_set_w_AverageBinPA_predictions = ...
    [unknown_with_features(:,1:27),...
    table(Final_Result_prediction,'VariableNames',{'AverageBinPA_Prediction'}),...
    unknown_with_features(:,28:end)];

writetable(Unknown_data_set_w_AverageBinPA_predictions, 'Unknown_data_set_w_AverageBinPA_predictions.xlsx','Sheet',1);

