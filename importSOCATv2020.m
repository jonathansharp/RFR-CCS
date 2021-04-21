disp('Importing SOCATv2020 observations');

%% Load CCS SOCATv2020 data
load('/Volumes/2TB Hard Drive/SOCAT/SOCATv2020_lat_15_60_lon_220_255.mat');
% Assemble into structure
SOCATv2020.all.expocode = Expocode;
SOCATv2020.all.latitude = latitude;
SOCATv2020.all.longitude = longitude;
SOCATv2020.all.temperature = SST;
SOCATv2020.all.salinity = sal;
SOCATv2020.all.day = day;
SOCATv2020.all.month = mon;
SOCATv2020.all.year = yr;
SOCATv2020.all.hour = hh;
SOCATv2020.all.minute = mm;
SOCATv2020.all.second = ss;
SOCATv2020.all.pressure = NCEP_SLP;
SOCATv2020.all.dist_to_land = dist_to_land;
SOCATv2020.all.fCO2 = fCO2rec;
SOCATv2020.all.fCO2_flag = fCO2rec_flag;
SOCATv2020.all.fCO2_src = fCO2rec_src;
SOCATv2020.all.flag = QC_Flag;
% Clean up
clear allexts ans coastlat coastlon day dotsize Expocode fCO2rec fCO2rec_flag
clear fCO2rec_src fname fpath GVCO2 hcb hdrline hdrs hh i latitude longitude
clear mess mm mon mx plotv plotvc projDir QC_Flag runcode sal ss SST StartText
clear x1lim xCO2water_equ_dry xx y1lim yr yt ytl yy hhh dist_to_land map_height

%% Convert fCO2 to pCO2
% Calculate fugacity factor
TempK = SOCATv2020.all.temperature+273.15;
Delta = (57.7 - 0.118.*TempK);
b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;
% For a mixture of CO2 and air at 1 atm (at low CO2 concentrations),
P1atm = SOCATv2020.all.pressure./1000; % in bar
RGasConstant = 83.14462618;
FugFac = exp((b + 2.*Delta).*P1atm./(RGasConstant.*TempK));
% Convert fCO2 to pCO2
SOCATv2020.all.pCO2 = SOCATv2020.all.fCO2./FugFac;

clear TempK Delta b P1atm RGasConstant FugFac

%% Remove observations older than 1998
idxyr = SOCATv2020.all.year >= 1998;
vars = fieldnames(SOCATv2020.all);
for n = 1:numel(vars)
        tempvar = SOCATv2020.all.(string(vars(n)));
        SOCATv2020.all.(string(vars(n))) = tempvar(idxyr);
end

clear idxyr vars

%% Remove observations with flags other than 2 and A/B/C/D
idxflag = SOCATv2020.all.fCO2_flag == 2 & ...
    (strcmp(SOCATv2020.all.flag,'A') | strcmp(SOCATv2020.all.flag,'B') | ...
    strcmp(SOCATv2020.all.flag,'C') | strcmp(SOCATv2020.all.flag,'D'));
vars = fieldnames(SOCATv2020.all);
for n = 1:numel(vars)
        tempvar = SOCATv2020.all.(string(vars(n)));
        SOCATv2020.all.(string(vars(n))) = tempvar(idxflag);
end

clear idxflag vars tempvar

%% Obtain unique integers for each expocode
SOCATv2020.all.cruise = nan(size(SOCATv2020.all.expocode));
cruiselist=unique(SOCATv2020.all.expocode);
for n=1:numel(cruiselist)
    idx = strcmp(cruiselist(n),SOCATv2020.all.expocode);
    SOCATv2020.all.cruise(idx) = n;
end

clear n cruiselist idx

%% Determine months since Jan 1 1998
SOCATv2020.all.month_since_1998 = (SOCATv2020.all.year-1998).*12 + SOCATv2020.all.month;

% Visualize number of observations
figure
histogram(SOCATv2020.all.month_since_1998);
ylabel('Number of observations');
xlabel('Months since January 1, 1998');

% %% Plot observations by date
% SOCATv2020.all.date = datenum([SOCATv2020.all.year SOCATv2020.all.month SOCATv2020.all.day]);
% 
% % Visualize number of observations
% figure
% histogram(SOCATv2020.all.date);
% ylabel('Number of observations');
% xlabel('Year');
% ylim([0 2e4]);
% datetick('x')

if split==1

%% Remove observations randomly based on platform
disp('Splitting into training (80%) and test (20%) data randomly');
for h = 1:numsplits
    if h == 1
        rng(rng_seed);
    end
    rng(randi([1 100]));
    plat = unique(SOCATv2020.all.expocode);
    numplats = numel(plat);
    pertest = 0.2;
    test_plats = randperm(numplats,round(pertest*numplats))';
    test_idx = ismember(SOCATv2020.all.expocode,plat(test_plats));
    vars = fieldnames(SOCATv2020.all);
    for n = 1:numel(vars)
            tempvar = SOCATv2020.all.(string(vars(n)));
            SOCATv2020.(strcat('train',num2str(h))).(string(vars(n))) = tempvar(~test_idx);
            SOCATv2020.(strcat('test',num2str(h))).(string(vars(n))) = tempvar(test_idx);
    end
end

clear n plat numplats pertest test_plats test_idx vars tempvar

elseif split==2

%% Remove observations based on certain moorings
disp('Splitting into training and test data by removing specific moorings (WA and CCE1)');
plat = unique(SOCATv2020.all.expocode);
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
test_idx = ismember(SOCATv2020.all.expocode,test_plat);
vars = fieldnames(SOCATv2020.all);
for n = 1:numel(vars)
        tempvar = SOCATv2020.all.(string(vars(n)));
        SOCATv2020.train1.(string(vars(n))) = tempvar(~test_idx);
        SOCATv2020.test1.(string(vars(n))) = tempvar(test_idx);
end

clear n tempvar plat moor_idx moor_type elim_idx test_plat test_idx vars

elseif split==3
    
%% Remove observations based on year
for h = 1:numel(startyear)
test_idx = ismember(SOCATv2020.all.year,startyear(h):5:2020);
vars = fieldnames(SOCATv2020.all);
for n = 1:numel(vars)
    tempvar = SOCATv2020.all.(string(vars(n)));
    SOCATv2020.(strcat('train',num2str(h))).(string(vars(n))) = tempvar(~test_idx);
    SOCATv2020.(strcat('test',num2str(h))).(string(vars(n))) = tempvar(test_idx);
end

clear n test_idx vars tempvar

end

end

clear a
