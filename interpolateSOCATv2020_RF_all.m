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

% Import mooring data
load('/Users/sharp/Documents/DATA/MOORINGS/mooring_data_processed.mat');
Y_moor = [];
idxlon = [];
idxlat = [];
idxmon = [];
moornames = fieldnames(MOORING);
for k=1:numel(moornames)
    % Determine mooring month since 1998
    MOORING.(moornames{k}).month_since_1998 = ...
        (MOORING.(moornames{k}).date(:,1)-1998).*12 + MOORING.(moornames{k}).date(:,2);
    % Determine means of each month since 1998
    for m=1:max(SOCATv2020_grid.month_since_1998)
        Y_moor = [Y_moor;mean(MOORING.(moornames{k}).pCO2SW(MOORING.(moornames{k}).month_since_1998==m),'omitnan')];
        idxlon = [idxlon;find(abs(SOCATv2020_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{k}).lon)) == ...
            min(min(abs(SOCATv2020_grid.longitude(:,1,1)-nanmean(MOORING.(moornames{k}).lon)))))];
        idxlat = [idxlat;find(abs(SOCATv2020_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{k}).lat)) == ...
            min(min(abs(SOCATv2020_grid.latitude(1,:,1)-nanmean(MOORING.(moornames{k}).lat)))))];
        idxmon = [idxmon;m];
    end
end
X_moor = [];
for n = 1:size(Y_moor,1)
X_moor = [X_moor;...
     SOCATv2020_grid.distance_from_shore(idxlon(n),idxlat(n),idxmon(n)) ...
     SOCATv2020_grid.SSS(idxlon(n),idxlat(n),idxmon(n)) ...
     SOCATv2020_grid.SSH(idxlon(n),idxlat(n),idxmon(n)) ...
     SOCATv2020_grid.SST(idxlon(n),idxlat(n),idxmon(n)) ...
     log10(SOCATv2020_grid.CHL(idxlon(n),idxlat(n),idxmon(n))) ...
     log10(SOCATv2020_grid.MLD(idxlon(n),idxlat(n),idxmon(n))) ...
     SOCATv2020_grid.wind_speed(idxlon(n),idxlat(n),idxmon(n)) ...
     SOCATv2020_grid.pCO2_atm(idxlon(n),idxlat(n),idxmon(n)) ...
     SOCATv2020_grid.year_since_1997(idxlon(n),idxlat(n),idxmon(n)) ...
     SOCATv2020_grid.month_sin(idxlon(n),idxlat(n),idxmon(n)) ...
     SOCATv2020_grid.month_cos(idxlon(n),idxlat(n),idxmon(n))];
end

%% Predictor variable combinations
Index_var = true(size(X,2)+1,size(X,2));
for h = 1:size(X,2)
Index_var(h+1,h) = false;
end

clear Table Table_row_headers
for v = 1:1%size(X,2)+1 % Run for different sets of predictor variables

headers_loop = headers(~Index_var(v,:));

%% ALL Data
disp('**Can uncomment this section to use all data to construct model')
% Index training data based on test data and available variables
Index_all = ~isnan(SOCATv2020_grid.all.pco2_ave_weighted) & ...
            ~isnan(SOCATv2020_grid.SSS) & ...
            ~isnan(SOCATv2020_grid.MLD) & ...
            ~isnan(SOCATv2020_grid.SST) & ...
            ~isnan(SOCATv2020_grid.SSH) & ...
            ~isnan(SOCATv2020_grid.CHL);
% Predictor variables:
X_train = X(Index_all,:);
X_train = X_train(:,Index_var(v,:));
X_test = X(Index_all,:);
X_test = X_test(:,Index_var(v,:));
X_moor_all = X_moor(:,Index_var(v,:));
% Response variable (pCO2):
Y_train = SOCATv2020_grid.all.pco2_ave_weighted(:);
Y_train = Y_train(Index_all);
Y_test = SOCATv2020_grid.all.pco2_ave_weighted(:);
Y_test = Y_test(Index_all);

%% Train Random Forest Model

% Set parameters for RF model
% rng(1); % For reproducibility
nTrees = 1000;
minLeafSize = 2;
numpredictors = 6;

% Leaf size test
% disp('**Can uncomment this section to try different values of minLeafSize')
% leaf = [1 2 3 4 5 10 20];
% col = 'rbgcmyk';
% figure
% for i=1:length(leaf)
%     b = TreeBagger(nTrees,X_train,Y_train,'Method','regression','OOBPrediction','on',...
%         'MinLeafSize',leaf(i),'NumPredictorsToSample',numpredictors);
%     plot(oobError(b),col(i))
%     hold on
% end
% xlabel('Number of Grown Trees')
% ylabel('Mean Squared Error') 
% legend({'1' '2' '3' '4' '5' '10' '20'},'Location','NorthEast')
% hold off

% Construct RF model
b = TreeBagger(nTrees,X_train,Y_train,'Method','regression','OOBPrediction','on',...
    'MinLeafSize',minLeafSize,'NumPredictorsToSample',numpredictors);

% Plot Out-of-Bag MSE based on tree number
figure
plot(oobError(b),'k','LineWidth',2);
xlabel('Number of Grown Trees');
ylabel('Out-of-Bag Mean Squared Error');

% Plot importance of each predictor
% figure
% bar(b.OOBPermutedPredictorDeltaError);
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

RF_pred_moor = predict(b,X_moor_all);
mdl_moor = fitlm(Y_moor,RF_pred_moor);

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

% Plot measured vs. predicted pCO2 values
figure; hold on;
set(gca,'fontsize',16);
colormap(cmocean('thermal'));
scatter_kde(Y_moor, RF_pred_moor, 'filled', 'MarkerSize', 15); 
plot([100 1000],[100 1000],'k--','LineWidth',2);
xlabel('Mooring measured monthly mean pCO2','fontsize',16);
xlim([100 1000]);
ylabel('RFR-estimated pCO_{2}','fontsize',16);
ylim([100 1000]);
c=colorbar; 
c.Label.String = 'Relative data density';

% Determine error statistics
err_RF = RF_pred-Y_test; % Errors
rmse_RF = sqrt(mean(err_RF.^2)); % RMSE
err_moor = RF_pred_moor-Y_moor; % Mooring Errors
rmse_moor = sqrt(mean(err_moor.^2,'omitnan')); % Mooring RMSE
out_of_bag_error = oobError(b);
fprintf('Mean of RF prediction: \n %f \n',nanmean(err_RF));
fprintf('RMSE from RF prediction: \n %f \n',rmse_RF);
fprintf('R^2 of SOCAT~Predicted: \n %f \n',mdl.Rsquared.Ordinary);
fprintf('RMSE from mooring predictions: \n %f \n',rmse_moor);
fprintf('OOB RMSE: \n %f \n',sqrt(out_of_bag_error(end)));

% Plot histogram
% edges = -200:0.5:200;
% figure; histogram(err_RF,edges); xlabel('Error'); ylabel('Number of Points');
% title('RF');

Table_row_headers{v,:} = strjoin(headers_loop);
Table(v,1) = nanmean(err_RF);
Table(v,2) = rmse_RF;
Table(v,3) = mdl.Rsquared.Ordinary;
Table(v,4) = rmse_moor;
Table(v,5) = sqrt(out_of_bag_error(end));

end

if grid == 1

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
SOCATv2020_grid.pco2_RF = nan(size(SOCATv2020_grid.latitude,1)*...
                                   size(SOCATv2020_grid.latitude,2)*...
                                   size(SOCATv2020_grid.latitude,3),1);
% Fill and reshape predicted pCO2 values to 3D matrix
SOCATv2020_grid.pco2_RF(Index_grid) = RF_pred_grid;
SOCATv2020_grid.pco2_RF = ...
    reshape(SOCATv2020_grid.pco2_RF,size(SOCATv2020_grid.latitude));

%% Remove trend from time series using monthly domain mean pCO2
% Calculate domain means
disp('Removing trend from pCO2');
SOCATv2020_grid.pco2_RF_dom_mean = squeeze(nanmean(nanmean(SOCATv2020_grid.pco2_RF,1),2));
SOCATv2020_grid.pco2_RF_dom_std = nanstd(squeeze(nanstd(SOCATv2020_grid.pco2_RF,1)),1)';
% Fit trends
[SOCATv2020_grid.pco2_RF_yf,~,~] = leastsq(1:max(SOCATv2020_grid.month_since_1998),SOCATv2020_grid.pco2_RF_dom_mean,0,0,0);
% Remove difference from mean for each month:
for m = 1:max(SOCATv2020_grid.month_since_1998)
    SOCATv2020_grid.pco2_RF_mon_mean_detrend(:,:,m) = SOCATv2020_grid.pco2_RF(:,:,m) + (nanmean(SOCATv2020_grid.pco2_RF_yf) - SOCATv2020_grid.pco2_RF_yf(m));
end
clear m yr x

%% Determine monthly means for one annual cycle (1998-2019 climatology)
disp('Determining monthly means for one annual cycle');
% Pre-allocate
SOCATv2020_grid.pco2_RF_clim = ...
    nan(size(SOCATv2020_grid.pco2_RF,1),...
    size(SOCATv2020_grid.pco2_RF,2),12);
SOCATv2020_grid.pco2_RF_clim_std = ...
    nan(size(SOCATv2020_grid.pco2_RF,1),...
    size(SOCATv2020_grid.pco2_RF,2),12);
SOCATv2020_grid.pco2_RF_clim_2015 = ...
    nan(size(SOCATv2020_grid.pco2_RF,1),...
    size(SOCATv2020_grid.pco2_RF,2),12);
for m = 1:12
    SOCATv2020_grid.pco2_RF_clim(:,:,m) = ...
        nanmean(SOCATv2020_grid.pco2_RF_mon_mean_detrend(:,:,m:12:end),3);
    SOCATv2020_grid.pco2_RF_clim_std(:,:,m) = ...
        std(SOCATv2020_grid.pco2_RF_mon_mean_detrend(:,:,m:12:end),0,3);
    SOCATv2020_grid.pco2_RF_clim_2015(:,:,m) = ...
        nanmean(SOCATv2020_grid.pco2_RF_mon_mean_detrend(:,:,m:12:end-48),3);
    SOCATv2020_grid.pco2_RF_clim_std_2015(:,:,m) = ...
        std(SOCATv2020_grid.pco2_RF_mon_mean_detrend(:,:,m:12:end-48),0,3);
end
clear m

%% Determine annual mean and seasonal amplitude
disp('Determining annual means and amplitudes');
% Determine annual mean values
SOCATv2020_grid.pco2_RF_annmean = nanmean(SOCATv2020_grid.pco2_RF_clim,3);
SOCATv2020_grid.pco2_RF_annmean_2015 = nanmean(SOCATv2020_grid.pco2_RF_clim_2015,3);
% Determine seasonal amplitudes
SOCATv2020_grid.pco2_RF_amp = ...
    max(SOCATv2020_grid.pco2_RF_clim,[],3,'omitnan') - ...
    min(SOCATv2020_grid.pco2_RF_clim,[],3,'omitnan');
SOCATv2020_grid.pco2_RF_amp_2015 = ...
    max(SOCATv2020_grid.pco2_RF_clim_2015,[],3,'omitnan') - ...
    min(SOCATv2020_grid.pco2_RF_clim_2015,[],3,'omitnan');
clear idx_month

end
