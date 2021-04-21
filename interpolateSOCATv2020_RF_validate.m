%% Random Forest estimation
disp('Train Random Forest model to estimate pCO2');

% Year relative to 1997:
SOCATv2020_grid.year_since_1997 = SOCATv2020_grid.year-1997;
% Transform month by sine and cosine:
SOCATv2020_grid.month_sin = sin((2.*pi.*SOCATv2020_grid.month)./12);
SOCATv2020_grid.month_cos = cos((2.*pi.*SOCATv2020_grid.month)./12);
% Remove bathymetry above sea level:
% SOCATv2020_grid.bottomdepth(SOCATv2020_grid.bottomdepth<0) = NaN;

%% Predictor variables
X = [SOCATv2020_grid.distance_from_shore(:) ...
     SOCATv2020_grid.SSS(:) ...
     SOCATv2020_grid.SSH(:) ...
     SOCATv2020_grid.SST(:) ...
     log10(SOCATv2020_grid.CHL(:)) ...
     log10(SOCATv2020_grid.MLD(:)) ...
     SOCATv2020_grid.wind_speed(:) ...
     SOCATv2020_grid.pCO2_atm(:) ...
     SOCATv2020_grid.year_since_1997(:) ...
     SOCATv2020_grid.month_sin(:) ...
     SOCATv2020_grid.month_cos(:)];
% Standard case:
% X = [SOCATv2020_grid.distance_from_shore(:) ...
%      SOCATv2020_grid.SSS(:) ...
%      SOCATv2020_grid.SSH(:) ...
%      SOCATv2020_grid.SST(:) ...
%      log10(SOCATv2020_grid.CHL(:)) ...
%      log10(SOCATv2020_grid.MLD(:)) ...
%      SOCATv2020_grid.wind_speed(:) ...
%      SOCATv2020_grid.pCO2_atm(:) ...
%      SOCATv2020_grid.year_since_1997(:) ...
%      SOCATv2020_grid.month_sin(:) ...
%      SOCATv2020_grid.month_cos(:) ...
%      SOCATv2020_grid.bottomdepth(:)];

% Predictor variable headers
headers = {'Dist' 'SSS' 'SSH' 'SST' 'CHL' 'MLD' 'Wind' 'aCO2' 'Year' 'Month'};
% Standard case:
% headers = {'Dist' 'SSS' 'SSH' 'SST' 'CHL' 'MLD' 'Wind' 'aCO2' 'Year' 'Month' 'Bathy'};

% %% Predictor variable combinations
% Index_var = true(size(X,2)+1,size(X,2));
% for h = 1:size(X,2)
% Index_var(h+1,h) = false;
% end

%Index_var = true(numel(startyear),size(X,2));
Index_var = true(numsplits,size(X,2));

%for v = 1:numsplits % Run for different splits
%for v = 1:numel(startyear) % Run for different splits
for v = 1:1 % Run once

headers_loop = headers(~Index_var(v,:));

%% Training Data
% Index training data based on training data and available variables
Index_train = ~isnan(SOCATv2020_grid.(strcat('train',num2str(v))).pco2_ave_weighted(:)) & ...
              ~isnan(SOCATv2020_grid.SSS(:)) & ...
              ~isnan(SOCATv2020_grid.SSH(:)) & ...
              ~isnan(SOCATv2020_grid.MLD(:)) & ...
              ~isnan(SOCATv2020_grid.CHL(:)) & ...
              ~isnan(SOCATv2020_grid.SST(:));
% Predictor variables:
X_train = X(Index_train,:);
X_train = X_train(:,Index_var(v,:));
% Response variable (pCO2):
Y_train = SOCATv2020_grid.(strcat('train',num2str(v))).pco2_ave_weighted(:);
Y_train = Y_train(Index_train);

%% Test Data
% Index training data based on training data and available variables
Index_test = ~isnan(SOCATv2020_grid.(strcat('test',num2str(v))).pco2_ave_weighted(:)) & ...
              ~isnan(SOCATv2020_grid.SSS(:)) & ...
              ~isnan(SOCATv2020_grid.SSH(:)) & ...
              ~isnan(SOCATv2020_grid.MLD(:)) & ...
              ~isnan(SOCATv2020_grid.CHL(:)) & ...
              ~isnan(SOCATv2020_grid.SST(:));
% Predictor variables:
X_test = X(Index_test,:);
X_test = X_test(:,Index_var(v,:));
% Response variable (pCO2):
Y_test = SOCATv2020_grid.(strcat('test',num2str(v))).pco2_ave_weighted(:);
Y_test = Y_test(Index_test);

%% Train Random Forest Model

% Set parameters for RF model
% rng(1); % For reproducibility
nTrees = 1000;
minLeafSize = 2;
numpredictors = 6;

% Leaf size test
disp('**Can uncomment this section to try different values of minLeafSize')
% leaf = [2 5 10 20 50];
% col = 'rbcmy';
% figure
% for i=1:length(leaf)
%     b = TreeBagger(nTrees,X,Y,'Method','R','OOBPrediction','On',...
% 			'MinLeafSize',leaf(i));
%     plot(oobError(b),col(i))
%     hold on
% end
% xlabel('Number of Grown Trees')
% ylabel('Mean Squared Error') 
% legend({'2' '5' '10' '20' '50'},'Location','NorthEast')
% hold off

% Construct RF model
b = TreeBagger(nTrees,X_train,Y_train,'Method','regression',...
    'MinLeafSize',minLeafSize,'NumPredictorsToSample',numpredictors);

% Plot Out-of-Bag MSE based on tree number
% figure
% plot(oobError(b),'k','LineWidth',2);
% xlabel('Number of Grown Trees');
% ylabel('Out-of-Bag Mean Squared Error');

% % Plot importance of each predictor
% figure
% bar(oobPermutedPredictorImportance(b));
% xlabel('Predictor') ;
% ylabel('Out-of-Bag Feature Importance');
% xticklabels(headers);

clear Index minLeafSize minLeafSize

%% Evaluate Random Forest model with test data
disp('Evaluate random forest model');

% Predict pCO2 for each test grid cell
RF_pred = predict(b,X_test);
% Fit linear relationship between gridded and predicted pCO2
mdl = fitlm(Y_test,RF_pred);

% Plot measured vs. predicted pCO2 values
figure; hold on;
set(gca,'fontsize',16);
colormap(cmocean('thermal'));
scatter_kde(Y_test, RF_pred, 'filled', 'MarkerSize', 15); 
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('Gridded pCO_{2} (test data)','fontsize',16);
xlim([100 1000]);
ylabel('RFR-estimated pCO_{2}','fontsize',16);
ylim([100 1000]);
c=colorbar; 
c.Label.String = 'Relative data density';

% Determine error statistics
err_RF = RF_pred-Y_test; % Errors
rmse_RF = sqrt(mean(err_RF.^2)); % RMSE
fprintf('Mean of RF prediction: \n %f \n',nanmean(err_RF));
fprintf('RMSE from RF prediction: \n %f \n',rmse_RF);
fprintf('R^2 of SOCAT~Predicted: \n %f \n',mdl.Rsquared.Ordinary);


% Plot histogram
% edges = -200:0.5:200;
% figure; histogram(err_RF,edges); xlabel('Error'); ylabel('Number of Points');
% title('RF');

Table_row_headers{v,:} = strjoin(headers_loop);
Table_val(v,1) = nanmean(err_RF);
Table_val(v,2) = rmse_RF;
Table_val(v,3) = mdl.Rsquared.Ordinary;

end

if grid_val ==1

%% Predict pCO2 values on SOCAT grid using RF
disp('Estimating pCO2 on grid using Random Forest model');

% Index training data based on available variables
Index_grid = ~isnan(SOCATv2020_grid.SSS) & ...
             ~isnan(SOCATv2020_grid.SSH) & ...
             ~isnan(SOCATv2020_grid.MLD) & ...
             ~isnan(SOCATv2020_grid.CHL) & ...
             ~isnan(SOCATv2020_grid.SST);
% Predictor variables:
X_grid = X(Index_grid,:);
% Predict pCO2 for each grid cell
RF_pred_grid = predict(b,X_grid);

% Pre-allocate predicted pCO2
SOCATv2020_grid.pco2_RF_validate = nan(size(SOCATv2020_grid.latitude,1)*...
                                   size(SOCATv2020_grid.latitude,2)*...
                                   size(SOCATv2020_grid.latitude,3),1);
% Fill and reshape predicted pCO2 values to 3D matrix
SOCATv2020_grid.pco2_RF_validate(Index_grid) = RF_pred_grid;
SOCATv2020_grid.pco2_RF_validate = ...
    reshape(SOCATv2020_grid.pco2_RF_validate,size(SOCATv2020_grid.latitude));

%% Remove trend from time series using monthly domain mean pCO2
% Calculate domain means
disp('Removing trend from pCO2');
SOCATv2020_grid.pco2_RF_validate_dom_mean = squeeze(nanmean(nanmean(SOCATv2020_grid.pco2_RF_validate,1),2));
SOCATv2020_grid.pco2_RF_validate_dom_std = nanstd(squeeze(nanstd(SOCATv2020_grid.pco2_RF_validate,1)),1)';
% Fit trends
[SOCATv2020_grid.pco2_RF_validate_yf,~,~] = leastsq(1:max(SOCATv2020_grid.month_since_1998),SOCATv2020_grid.pco2_RF_validate_dom_mean,0,0,0);
% Remove difference from mean for each month:
for m = 1:max(SOCATv2020_grid.month_since_1998)
    SOCATv2020_grid.pco2_RF_validate_mon_mean_detrend(:,:,m) = SOCATv2020_grid.pco2_RF_validate(:,:,m) + (nanmean(SOCATv2020_grid.pco2_RF_validate_yf) - SOCATv2020_grid.pco2_RF_validate_yf(m));
end
clear m yr x

%% Determine monthly means for one annual cycle (1998-2020 climatology)
disp('Determining monthly means for one annual cycle');
% Pre-allocate
SOCATv2020_grid.pco2_RF_validate_clim = ...
    nan(size(SOCATv2020_grid.pco2_RF_validate,1),...
    size(SOCATv2020_grid.pco2_RF_validate,2),12);
SOCATv2020_grid.pco2_RF_validate_clim_std = ...
    nan(size(SOCATv2020_grid.pco2_RF_validate,1),...
    size(SOCATv2020_grid.pco2_RF_validate,2),12);
for m = 1:12
    SOCATv2020_grid.pco2_RF_validate_clim(:,:,m) = ...
        nanmean(SOCATv2020_grid.pco2_RF_validate_mon_mean_detrend(:,:,m:12:end),3);
    SOCATv2020_grid.pco2_RF_validate_clim_std(:,:,m) = ...
        std(SOCATv2020_grid.pco2_RF_validate_mon_mean_detrend(:,:,m:12:end),0,3);
end
clear m

%% Determine annual mean and seasonal amplitude
disp('Determining annual means and amplitudes');
% Determine annual mean values
SOCATv2020_grid.pco2_RF_validate_annmean = nanmean(SOCATv2020_grid.pco2_RF_validate_clim,3);
% Determine seasonal amplitudes
SOCATv2020_grid.pco2_RF_validate_amp = ...
    max(SOCATv2020_grid.pco2_RF_validate_clim,[],3,'omitnan') - ...
    min(SOCATv2020_grid.pco2_RF_validate_clim,[],3,'omitnan');
clear idx_month

end
