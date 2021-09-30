disp('Importing SOCATv2021 observations');

%% Load CCS SOCATv2021 data
load('Data/SOCATv2021_lat_15_60_lon_220_255.mat');
% Assemble into structure
SOCATv2021.all.expocode = Expocode;
SOCATv2021.all.latitude = latitude;
SOCATv2021.all.longitude = longitude;
SOCATv2021.all.temperature = SST;
SOCATv2021.all.salinity = sal;
SOCATv2021.all.day = day;
SOCATv2021.all.month = mon;
SOCATv2021.all.year = yr;
SOCATv2021.all.hour = hh;
SOCATv2021.all.minute = mm;
SOCATv2021.all.second = ss;
SOCATv2021.all.pressure = NCEP_SLP;
SOCATv2021.all.dist_to_land = dist_to_land;
SOCATv2021.all.fCO2 = fCO2rec;
SOCATv2021.all.fCO2_flag = fCO2rec_flag;
SOCATv2021.all.fCO2_src = fCO2rec_src;
SOCATv2021.all.flag = QC_Flag;
% Clean up
clear allexts ans coastlat coastlon day dist_to_land dotsize ETOPO2_depth
clear Expocode fCO2rec fCO2rec_flag fCO2rec_src fCO2water_equ_wet fCO2water_SST_wet
clear fname fpath GVCO2 hcb hdrline hdrs hh hhh i latitude longitude map_height
clear mess mm mon mx NCEP_SLP pCO2water_equ_wet pCO2water_SST_wet Pequ plotv
clear plotvc PPPP projDir QC_Flag runcode sal sample_depth Source_DOI ss
clear SST StartText Tequ version WOA_SSS x1lim xCO2water_equ_dry xCO2water_SST_dry
clear xx y1lim yr yt ytl yy

%% Convert fCO2 to pCO2
% Calculate fugacity factor
TempK = SOCATv2021.all.temperature+273.15;
Delta = (57.7 - 0.118.*TempK);
b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;
P1atm = SOCATv2021.all.pressure./1000; % in bar
RGasConstant = 83.14462618;
FugFac = exp((b + 2.*Delta).*P1atm./(RGasConstant.*TempK));
% Convert fCO2 to pCO2
SOCATv2021.all.pCO2 = SOCATv2021.all.fCO2./FugFac;

%Vapor Pressure Correction????
%SOCATv2021.all.vapor_pressure = vpress(SOCATv2021.all.salinity,SOCATv2021.all.temperature);

clear TempK Delta b P1atm RGasConstant FugFac

%% Remove observations older than 1998
idxyr = SOCATv2021.all.year >= 1998;
vars = fieldnames(SOCATv2021.all);
for n = 1:numel(vars)
        tempvar = SOCATv2021.all.(string(vars(n)));
        SOCATv2021.all.(string(vars(n))) = tempvar(idxyr);
end

clear n tempvar idxyr vars

%% Remove observations with flags other than 2 and A/B/C/D
idxflag = SOCATv2021.all.fCO2_flag == 2 & ...
    (strcmp(SOCATv2021.all.flag,'A') | strcmp(SOCATv2021.all.flag,'B') | ...
    strcmp(SOCATv2021.all.flag,'C') | strcmp(SOCATv2021.all.flag,'D'));
vars = fieldnames(SOCATv2021.all);
for n = 1:numel(vars)
        tempvar = SOCATv2021.all.(string(vars(n)));
        SOCATv2021.all.(string(vars(n))) = tempvar(idxflag);
end

clear n idxflag vars tempvar

%% Obtain unique integers for each expocode
SOCATv2021.all.cruise = nan(size(SOCATv2021.all.expocode));
cruiselist=unique(SOCATv2021.all.expocode);
for n=1:numel(cruiselist)
    idx = strcmp(cruiselist(n),SOCATv2021.all.expocode);
    SOCATv2021.all.cruise(idx) = n;
end

clear n cruiselist idx

%% Determine months since Jan 1 1998
SOCATv2021.all.month_since_1998 = (SOCATv2021.all.year-1998).*12 + SOCATv2021.all.month;

% %% Visualize number of observations
% figure
% histogram(SOCATv2021.all.month_since_1998);
% ylabel('Number of observations');
% xlabel('Months since January 1, 1998');

% %% Plot observations by date
% SOCATv2021.all.date = datenum([SOCATv2021.all.year SOCATv2021.all.month SOCATv2021.all.day]);
% 
% % Visualize number of observations
% figure
% histogram(SOCATv2021.all.date);
% ylabel('Number of observations');
% xlabel('Year');
% ylim([0 2e4]);
% datetick('x')


%% Remove observations randomly based on platform
if split==1
disp('Splitting into training (80%) and test (20%) data randomly');
for h = 1:numsplits
    plat = unique(SOCATv2021.all.expocode);
    numplats = numel(plat);
    pertest = 0.2;
    test_plats = randperm(numplats,round(pertest*numplats))';
    test_idx = ismember(SOCATv2021.all.expocode,plat(test_plats));
    vars = fieldnames(SOCATv2021.all);
    for n = 1:numel(vars)
            tempvar = SOCATv2021.all.(string(vars(n)));
            SOCATv2021.(strcat('train',num2str(h))).(string(vars(n))) = tempvar(~test_idx);
            SOCATv2021.(strcat('test',num2str(h))).(string(vars(n))) = tempvar(test_idx);
    end
end

%SOCATv2021 = rmfield(SOCATv2021,'all');

clear n plat numplats pertest test_plats test_idx vars tempvar

%% Remove observations based on certain moorings
elseif split==2
disp('Splitting into training and test data by removing all moorings');
omitmoors = {'WA' 'CCE1' 'CCE2' 'SEAK' 'NH10' 'CB06' 'LaPush' 'KwakshuaChannel' 'Dabob' 'Exp'};
plat = unique(SOCATv2021.all.expocode);
moor_idx = ~cellfun(@isempty,regexp(plat,'3164')); % This determines which expocodes are moorings
moor_idx = moor_idx | ~cellfun(@isempty,regexp(plat,'187F'));
moor_type = {'Kwakshua';'WA';'WA';'WA';'CCE1';'CCE1';'WA';'CCE1';'CCE2';'WA';'LaPush';...
             'CCE1';'CCE2';'CCE1';'LaPush';'CCE2';'LaPush';'WA';'SEAK';'CCE2';'LaPush';...
             'Exp';'WA';'Exp';'CCE1';'LaPush';'SEAK';'Exp';'NH10';'CCE2';'Exp';...
             'SEAK';'LaPush';'WA'; 'Exp';'CCE1';'SEAK';'CCE2';'LaPush';'WA';'NH10';...
             'CCE1';'Dabob';'LaPush'; 'CCE2';'NH10';'CCE1';'LaPush'; 'CCE2';'LaPush';'CB06';...
             'WA';'WA';'CB06';'CB06'}; % This represents the identity of each mooring
elim_idx = ismember(moor_type,omitmoors); % This chooses moorings to eliminate
test_plat = plat(moor_idx);
test_plat = test_plat(elim_idx);
test_idx = ismember(SOCATv2021.all.expocode,test_plat);
vars = fieldnames(SOCATv2021.all);
for n = 1:numel(vars)
        tempvar = SOCATv2021.all.(string(vars(n)));
        SOCATv2021.train1.(string(vars(n))) = tempvar(~test_idx);
        SOCATv2021.test1.(string(vars(n))) = tempvar(test_idx);
end

%SOCATv2021 = rmfield(SOCATv2021,'all');

clear n tempvar plat moor_idx moor_type elim_idx test_plat test_idx vars
    
%% Remove observations based on year
elseif split==3
for h = 1:numel(startyear)
test_idx = ismember(SOCATv2021.all.year,startyear(h):5:2020);
vars = fieldnames(SOCATv2021.all);
for n = 1:numel(vars)
    tempvar = SOCATv2021.all.(string(vars(n)));
    SOCATv2021.(strcat('train',num2str(h))).(string(vars(n))) = tempvar(~test_idx);
    SOCATv2021.(strcat('test',num2str(h))).(string(vars(n))) = tempvar(test_idx);
end

%SOCATv2021 = rmfield(SOCATv2021,'all');

clear n test_idx vars tempvar

end

end

clear a
