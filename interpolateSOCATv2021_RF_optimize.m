%% Random Forest estimation
disp('Train Random Forest model to estimate pCO2');
clear Table_row_headers Table_val

% Year relative to 1997:
SOCATv2021_grid.year_since_1997 = SOCATv2021_grid.year-1997;
% Transform month by sine and cosine:
SOCATv2021_grid.month_sin = sin((2.*pi.*SOCATv2021_grid.month)./12);
SOCATv2021_grid.month_cos = cos((2.*pi.*SOCATv2021_grid.month)./12);
% Remove bathymetry above sea level:
% SOCATv2021_grid.bottomdepth(SOCATv2021_grid.bottomdepth<0) = NaN;

%% Predictor variables
X = [SOCATv2021_grid.distance_from_shore(:) ...
     SOCATv2021_grid.SSS(:) ...
     SOCATv2021_grid.SST(:) ...
     log10(SOCATv2021_grid.CHL(:)) ...
     log10(SOCATv2021_grid.MLD(:)) ...
     SOCATv2021_grid.wind_speed(:) ...
     SOCATv2021_grid.pCO2_atm(:) ...
     SOCATv2021_grid.year_since_1997(:) ...
     SOCATv2021_grid.month_sin(:) ...
     SOCATv2021_grid.month_cos(:)];
% Standard case:
% X = [SOCATv2021_grid.distance_from_shore(:) ...
%      SOCATv2021_grid.latitude(:) ...
%      SOCATv2021_grid.longitude(:) ...
%      SOCATv2021_grid.SSS(:) ...
%      SOCATv2021_grid.SSH(:) ...
%      SOCATv2021_grid.SST(:) ...
%      log10(SOCATv2021_grid.CHL(:)) ...
%      log10(SOCATv2021_grid.MLD(:)) ...
%      SOCATv2021_grid.wind_speed(:) ...
%      SOCATv2021_grid.pCO2_atm(:) ...
%      SOCATv2021_grid.year_since_1997(:) ...
%      SOCATv2021_grid.month_sin(:) ...
%      SOCATv2021_grid.month_cos(:) ...
%      SOCATv2021_grid.bottomdepth(:)];

% Predictor variable headers
headers = {'Dist' 'SSS' 'SST' 'CHL' 'MLD' 'Wind' 'aCO2' 'Year' 'sin(Mo.)' 'cos(Mo.)'};
% Standard case:
% headers = {'Dist' 'Lat.' 'Lon.' 'SSS' 'SSH' 'SST' 'CHL' 'MLD' 'Wind' 'aCO2' 'Year' 'sin(Mo.)' 'cos(Mo.)' 'Bathy'};

% %% Predictor variable combinations
Index_var = true(size(X,2)+1,size(X,2));
for h = 1:size(X,2)
Index_var(h+1,h) = false;
end

%% Split dataset into training, validation, and evaluation data
rng(30); % for reproducibility
val_per = 0.1; % percentage to use as validation data
eval_per = 0.1; % percentage to use as evaluation data
num_idx = randperm(length(X))';
val_idx = num_idx <= val_per.*length(X);
eval_idx = num_idx >  val_per.*length(X) & ...
           num_idx <= (val_per+eval_per).*length(X);
train_idx = num_idx >  (val_per+eval_per).*length(X);

% Data distribution by time
figure; hold on;
edges=datenum([[1998:2020]' ones(23,1) ones(23,1)]);
histogram(datenum([SOCATv2021_grid.year(train_idx) SOCATv2021_grid.month(train_idx) ones(sum(train_idx),1)]),edges);
histogram(datenum([SOCATv2021_grid.year(val_idx) SOCATv2021_grid.month(val_idx) ones(sum(val_idx),1)]),edges);
histogram(datenum([SOCATv2021_grid.year(eval_idx) SOCATv2021_grid.month(eval_idx) ones(sum(val_idx),1)]),edges);
legend({'Training' 'Validation' 'Evaluation'})
datetick('x'); xlim([datenum([1997 1 1]) datenum([2021 1 1])]);
ylabel('Data points within each year');

% Data distribution by latitude
figure; hold on;
edges=[15:5:60];
histogram(SOCATv2021_grid.latitude(train_idx),edges);
histogram(SOCATv2021_grid.latitude(val_idx),edges);
histogram(SOCATv2021_grid.latitude(eval_idx),edges);
legend({'Training' 'Validation' 'Evaluation'})
ylabel('Data points within each latitude range');

% Data distribution by longitude
figure; hold on;
edges=[210:5:255];
histogram(SOCATv2021_grid.longitude(train_idx),edges);
histogram(SOCATv2021_grid.longitude(val_idx),edges);
histogram(SOCATv2021_grid.longitude(eval_idx),edges);
legend({'Training' 'Validation' 'Evaluation'})
ylabel('Data points within each longitude range');

%% Train model
w = 1;
for v = 1%:size(Index_var,1)
% v = 1;
% for w = 1:6

%headers_loop = headers(~Index_var(v,:));
headers_loop = w;

%% Training Data
% Index training data based on training data and available variables
Index_train = ~isnan(SOCATv2021_grid.all.pco2_ave_weighted(:)) & ...
              ~isnan(SOCATv2021_grid.SSS(:)) & ...
              ~isnan(SOCATv2021_grid.SSH(:)) & ...
              ~isnan(SOCATv2021_grid.MLD(:)) & ...
              ~isnan(SOCATv2021_grid.CHL(:)) & ...
              ~isnan(SOCATv2021_grid.SST(:)) & ...
              train_idx;
% Predictor variables:
X_train = X(Index_train,:);
X_train = X_train(:,Index_var(v,:));
% Response variable (pCO2):
Y_train = SOCATv2021_grid.all.pco2_ave_weighted(:);
Y_train = Y_train(Index_train);

%% Validation Data
% Index training data based on training data and available variables
Index_val = ~isnan(SOCATv2021_grid.all.pco2_ave_weighted(:)) & ...
              ~isnan(SOCATv2021_grid.SSS(:)) & ...
              ~isnan(SOCATv2021_grid.SSH(:)) & ...
              ~isnan(SOCATv2021_grid.MLD(:)) & ...
              ~isnan(SOCATv2021_grid.CHL(:)) & ...
              ~isnan(SOCATv2021_grid.SST(:)) & ...
              val_idx;
% Predictor variables:
X_val = X(Index_val,:);
X_val = X_val(:,Index_var(v,:));
% Response variable (pCO2):
Y_val = SOCATv2021_grid.all.pco2_ave_weighted(:);
Y_val = Y_val(Index_val);

%% Evaluation Data
% Index training data based on training data and available variables
Index_eval = ~isnan(SOCATv2021_grid.all.pco2_ave_weighted(:)) & ...
              ~isnan(SOCATv2021_grid.SSS(:)) & ...
              ~isnan(SOCATv2021_grid.SSH(:)) & ...
              ~isnan(SOCATv2021_grid.MLD(:)) & ...
              ~isnan(SOCATv2021_grid.CHL(:)) & ...
              ~isnan(SOCATv2021_grid.SST(:)) & ...
              eval_idx;
% Predictor variables:
X_eval = X(Index_eval,:);
X_eval = X_eval(:,Index_var(v,:));
% Response variable (pCO2):
Y_eval = SOCATv2021_grid.all.pco2_ave_weighted(:);
Y_eval = Y_eval(Index_eval);

%% Train Random Forest Model

% Set parameters for RF model
rng(30); % For reproducibility
nTrees = 1200;
minLeafSize = [1 2 5 10 20];
numpredictors = 1:9;
infrac = [0.5 0.6 0.7 0.8 0.9 1.0];

% Leaf size test
% disp('**Can uncomment this section to try different values of minLeafSize')
% leaf = [1 2 5 10 20];
% col = 'rbcmy';
% figure
% for i=1:length(leaf)
%     b = TreeBagger(nTrees,X_train,Y_train,'Method','regression','OOBPrediction','on',...
% 			'MinLeafSize',leaf(i));
%     plot(oobError(b),col(i))
%     hold on
% end
% xlabel('Number of Grown Trees')
% ylabel('Mean Squared Error') 
% legend({'2' '5' '10' '20' '50'},'Location','NorthEast')
% hold off

% Construct RF model
b = TreeBagger(nTrees,X_train,Y_train,'Method','regression','OOBPrediction','on','InBagFraction',infrac(6),...
    'MinLeafSize',minLeafSize(3),'NumPredictorsToSample',numpredictors(6),'OOBPredictorImportance','on');
% b = TreeBagger(nTrees,X_train,Y_train,'Method','regression','OOBPrediction','on',...
%     'MinLeafSize',minLeafSize,'OOBPredictorImportance','on');

% Plot Out-of-Bag MSE based on tree number
% figure
% plot(sqrt(oobError(b)),'k','LineWidth',2);
% xlabel('Number of Grown Trees');
% ylabel('Out-of-Bag Root Mean Squared Error');

% Plot importance of each predictor
if v == 1
    figure
    bar(b.OOBPermutedPredictorDeltaError);
    xlabel('Predictor') ;
    ylabel('Out-of-Bag Feature Importance');
    xticklabels(headers(Index_var(v,:)));
end

Table_val(w,5) = min(sqrt(oobError(b)),[],1);

clear Index minLeafSize minLeafSize

%% Evaluate Random Forest model with validation data
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Evaluate random forest model with validation data');

% Predict pCO2 for each test grid cell
RF_pred = predict(b,X_val);
% Fit linear relationship between gridded and predicted pCO2
mdl = fitlm(Y_val,RF_pred);

% Plot measured vs. predicted pCO2 values
% figure; hold on;
% set(gca,'fontsize',16);
% colormap(cmocean('thermal'));
% scatter_kde(Y_val, RF_pred, 'filled', 'MarkerSize', 15);
% plot([100 1000],[100 1000],'k--','LineWidth',2);
% xlabel('Gridded pCO_{2} (validation data)','fontsize',16);
% xlim([100 1000]);
% ylabel('RFR-estimated pCO_{2}','fontsize',16);
% ylim([100 1000]);
% c=colorbar;
% c.Label.String = 'Relative data density';

% Determine error statistics
Coastal_idx = ~isnan(RF_pred) & SOCATv2021_grid.distance_from_shore(Index_val) <= 400;
Open_idx = ~isnan(RF_pred) & SOCATv2021_grid.distance_from_shore(Index_val) > 400;
err_RF = RF_pred-Y_val; % Errors
err_RF_open = RF_pred(Open_idx)-Y_val(Open_idx); % Errors
err_RF_coast = RF_pred(Coastal_idx)-Y_val(Coastal_idx); % Errors
rmse_RF = sqrt(mean(err_RF.^2)); % RMSE
fprintf('Mean of RF prediction: \n %f \n',nanmean(err_RF));
fprintf('RMSE from RF prediction: \n %f \n',rmse_RF);
fprintf('R^2 of SOCAT~Predicted: \n %f \n',mdl.Rsquared.Ordinary);
fprintf('Mapping uncert. (open): \n %f \n',std(err_RF_open));
fprintf('Mapping uncert. (coast): \n %f \n',std(err_RF_coast));

% Plot histogram
% edges = -200:0.5:200;
% figure; histogram(err_RF,edges); xlabel('Error'); ylabel('Number of Points');
% title('RF');

% for variables test
% Table_row_headers{v,:} = strjoin(headers_loop);
% Table_val(v,1) = nanmean(err_RF);
% Table_val(v,2) = rmse_RF;
% Table_val(v,3) = mdl.Rsquared.Ordinary;
% Table_val(v,4) = median(err_RF);

% for number of predictors test
Table_row_headers(w) = headers_loop;
Table_val(w,1) = nanmean(err_RF);
Table_val(w,2) = rmse_RF;
Table_val(w,3) = mdl.Rsquared.Ordinary;
Table_val(w,4) = median(err_RF);

end

%% Evaluate Random Forest model with evaluation data
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Evaluate random forest model with evaluation data');

% Predict pCO2 for each test grid cell
RF_pred = predict(b,X_eval);
% Fit linear relationship between gridded and predicted pCO2
mdl = fitlm(Y_eval,RF_pred);

% Plot measured vs. predicted pCO2 values
figure; hold on;
set(gca,'fontsize',16);
colormap(cmocean('thermal'));
scatter_kde(Y_eval, RF_pred, 'filled', 'MarkerSize', 15);
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('Gridded pCO_{2} (evaluation data)','fontsize',16);
xlim([100 1000]);
ylabel('RFR-estimated pCO_{2}','fontsize',16);
ylim([100 1000]);
c=colorbar;
c.Label.String = 'Relative data density';

% Determine error statistics
Coastal_idx = ~isnan(RF_pred) & SOCATv2021_grid.distance_from_shore(Index_eval) <= 400;
Open_idx = ~isnan(RF_pred) & SOCATv2021_grid.distance_from_shore(Index_eval) > 400;
err_RF = RF_pred-Y_eval; % Errors
err_RF_open = RF_pred(Open_idx)-Y_eval(Open_idx); % Errors
err_RF_coast = RF_pred(Coastal_idx)-Y_eval(Coastal_idx); % Errors
rmse_RF = sqrt(mean(err_RF.^2)); % RMSE
fprintf('Mean of RF prediction: \n %f \n',nanmean(err_RF));
fprintf('RMSE from RF prediction: \n %f \n',rmse_RF);
fprintf('R^2 of SOCAT~Predicted: \n %f \n',mdl.Rsquared.Ordinary);
fprintf('Mapping uncert. (open): \n %f \n',std(err_RF_open));
fprintf('Mapping uncert. (coast): \n %f \n',std(err_RF_coast));
