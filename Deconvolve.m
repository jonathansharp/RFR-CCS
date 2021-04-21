%% Calculate carbonate system parameters from TA and pCO2 and do seasonal decomposition of pCO2 and carbonate system parameters

%% Calculate monthly and annual means of SSS and SST
SOCATv2020_grid.SSS_clim = nan(size(SOCATv2020_grid.pco2_RF_clim));
SOCATv2020_grid.SST_clim = nan(size(SOCATv2020_grid.pco2_RF_clim));
for n = 1:12
    SOCATv2020_grid.SSS_clim(:,:,n) = mean(SOCATv2020_grid.SSS(:,:,n:12:end),3);
    SOCATv2020_grid.SST_clim(:,:,n) = mean(SOCATv2020_grid.SST(:,:,n:12:end),3);
end
SOCATv2020_grid.SSS_annmean = mean(SOCATv2020_grid.SSS_clim,3);
SOCATv2020_grid.SST_annmean = mean(SOCATv2020_grid.SST_clim,3);

%% Calculate monthly and annual means of TA
depth = repmat(5,size(SOCATv2020_grid.latitude(:,:,1:12)));
lat = SOCATv2020_grid.latitude(:,:,1:12);
lon = SOCATv2020_grid.longitude(:,:,1:12);
[TA,TA_UNCER] = ...
    ESPER_LIR(1,[lon(:) lat(:) depth(:)],...
    [SOCATv2020_grid.SSS_clim(:) SOCATv2020_grid.SST_clim(:)],[1 2], ...
            'Equations',8);
SOCATv2020_grid.TA_ESPER_LIR_clim = reshape(TA.TA,size(lat));
SOCATv2020_grid.TA_ESPER_LIR_annmean = mean(SOCATv2020_grid.TA_ESPER_LIR_clim,3);

%% Calculate monthly and annual means of DIC
CARB = CO2SYS(SOCATv2020_grid.pco2_RF_clim(:),SOCATv2020_grid.TA_ESPER_LIR_clim(:),4,1,SOCATv2020_grid.SSS_clim(:),...
       SOCATv2020_grid.SST_clim(:),NaN,0,NaN,0,0,0,0,1,10,1,2,2); % pCO2-TA calc
CARB(CARB==-999)=NaN;
SOCATv2020_grid.DIC_clim = reshape(CARB(:,2),size(SOCATv2020_grid.pco2_RF_clim));
SOCATv2020_grid.pH_clim = reshape(CARB(:,3),size(SOCATv2020_grid.pco2_RF_clim));
SOCATv2020_grid.OmA_clim = reshape(CARB(:,18),size(SOCATv2020_grid.pco2_RF_clim));
SOCATv2020_grid.DIC_annmean = mean(SOCATv2020_grid.DIC_clim,3);

%% Determine annual means of carbonate system params for each grid cell
CARB = CO2SYS(SOCATv2020_grid.TA_ESPER_LIR_annmean(:),SOCATv2020_grid.DIC_annmean(:),1,2,SOCATv2020_grid.SSS_annmean(:),...
       SOCATv2020_grid.SST_annmean(:),NaN,0,NaN,0,0,0,0,1,10,1,2,2); % TA-DIC calc
CARB(CARB==-999)=NaN;
SOCATv2020_grid.pCO2_annmean = reshape(CARB(:,4),size(SOCATv2020_grid.pco2_RF_annmean));
SOCATv2020_grid.pH_annmean = reshape(CARB(:,3),size(SOCATv2020_grid.pco2_RF_annmean));
SOCATv2020_grid.OmA_annmean = reshape(CARB(:,18),size(SOCATv2020_grid.pco2_RF_annmean));

%% Determine sensitivities of pCO2, pH, and Omega for each grid cell (dX/dY)
% ta
SENS = derivnum('par1',SOCATv2020_grid.TA_ESPER_LIR_clim(:),SOCATv2020_grid.DIC_clim(:),1,2,SOCATv2020_grid.SSS_clim(:),...
       SOCATv2020_grid.SST_clim(:),NaN,0,NaN,0,0,0,0,1,10,1,2,2);
dpco2_ta = reshape(SENS(:,4),size(SOCATv2020_grid.pco2_RF_clim));
dH_ta = reshape(SENS(:,3),size(SOCATv2020_grid.pco2_RF_clim));
dOmA_ta = reshape(SENS(:,11),size(SOCATv2020_grid.pco2_RF_clim));
% dic
SENS = derivnum('par2',SOCATv2020_grid.TA_ESPER_LIR_clim(:),SOCATv2020_grid.DIC_clim(:),1,2,SOCATv2020_grid.SSS_clim(:),...
       SOCATv2020_grid.SST_clim(:),NaN,0,NaN,0,0,0,0,1,10,1,2,2);
dpco2_dic = reshape(SENS(:,4),size(SOCATv2020_grid.pco2_RF_clim));
dH_dic = reshape(SENS(:,3),size(SOCATv2020_grid.pco2_RF_clim));
dOmA_dic = reshape(SENS(:,11),size(SOCATv2020_grid.pco2_RF_clim));
% tmp
SENS = derivnum('t',SOCATv2020_grid.TA_ESPER_LIR_clim(:),SOCATv2020_grid.DIC_clim(:),1,2,SOCATv2020_grid.SSS_clim(:),...
       SOCATv2020_grid.SST_clim(:),NaN,0,NaN,0,0,0,0,1,10,1,2,2);
dpco2_tmp = reshape(SENS(:,4),size(SOCATv2020_grid.pco2_RF_clim));
dH_tmp = reshape(SENS(:,3),size(SOCATv2020_grid.pco2_RF_clim));
dOmA_tmp = reshape(SENS(:,11),size(SOCATv2020_grid.pco2_RF_clim));
% sal
SENS = derivnum('s',SOCATv2020_grid.TA_ESPER_LIR_clim(:),SOCATv2020_grid.DIC_clim(:),1,2,SOCATv2020_grid.SSS_clim(:),...
       SOCATv2020_grid.SST_clim(:),NaN,0,NaN,0,0,0,0,1,10,1,2,2);
dpco2_sal = reshape(SENS(:,4),size(SOCATv2020_grid.pco2_RF_clim));
dH_sal = reshape(SENS(:,3),size(SOCATv2020_grid.pco2_RF_clim));
dOmA_sal = reshape(SENS(:,11),size(SOCATv2020_grid.pco2_RF_clim));

%% Calculate pCO2 components using sensitivities
SOCATv2020_grid.pCO2_TA = ...
    SOCATv2020_grid.pCO2_annmean + (SOCATv2020_grid.TA_ESPER_LIR_clim-SOCATv2020_grid.TA_ESPER_LIR_annmean).*dpco2_ta;
SOCATv2020_grid.pCO2_DIC = ...
    SOCATv2020_grid.pCO2_annmean + (SOCATv2020_grid.DIC_clim-SOCATv2020_grid.DIC_annmean).*dpco2_dic;
SOCATv2020_grid.pCO2_TMP = ...
    SOCATv2020_grid.pCO2_annmean + (SOCATv2020_grid.SST_clim-SOCATv2020_grid.SST_annmean).*dpco2_tmp;
SOCATv2020_grid.pCO2_SAL = ...
    SOCATv2020_grid.pCO2_annmean + (SOCATv2020_grid.SSS_clim-SOCATv2020_grid.SSS_annmean).*dpco2_sal;

%% Calculate pH components using sensitivities
SOCATv2020_grid.pH_TA = ...
    -log10( 10.^-(SOCATv2020_grid.pH_annmean) + (SOCATv2020_grid.TA_ESPER_LIR_clim-SOCATv2020_grid.TA_ESPER_LIR_annmean).*dH_ta.*1e-9 );
SOCATv2020_grid.pH_DIC = ...
    -log10( 10.^-(SOCATv2020_grid.pH_annmean) + (SOCATv2020_grid.DIC_clim-SOCATv2020_grid.DIC_annmean).*dH_dic.*1e-9 );
SOCATv2020_grid.pH_TMP = ...
    -log10( 10.^-(SOCATv2020_grid.pH_annmean) + (SOCATv2020_grid.SST_clim-SOCATv2020_grid.SST_annmean).*dH_tmp.*1e-9 );
SOCATv2020_grid.pH_SAL = ...
    -log10( 10.^-(SOCATv2020_grid.pH_annmean) + (SOCATv2020_grid.SSS_clim-SOCATv2020_grid.SSS_annmean).*dH_sal.*1e-9 );

%% Calculate Omega components using sensitivities
SOCATv2020_grid.OmA_TA = ...
    SOCATv2020_grid.OmA_annmean + (SOCATv2020_grid.TA_ESPER_LIR_clim-SOCATv2020_grid.TA_ESPER_LIR_annmean).*dOmA_ta;
SOCATv2020_grid.OmA_DIC = ...
    SOCATv2020_grid.OmA_annmean + (SOCATv2020_grid.DIC_clim-SOCATv2020_grid.DIC_annmean).*dOmA_dic;
SOCATv2020_grid.OmA_TMP = ...
    SOCATv2020_grid.OmA_annmean + (SOCATv2020_grid.SST_clim-SOCATv2020_grid.SST_annmean).*dOmA_tmp;
SOCATv2020_grid.OmA_SAL = ...
    SOCATv2020_grid.OmA_annmean + (SOCATv2020_grid.SSS_clim-SOCATv2020_grid.SSS_annmean).*dOmA_sal;


 %% Create plots
% Lat-lon index (1)
latitude_check = 28;
longitude_check = -120;
latidx = find(min(abs(SOCATv2020_grid.latitude(1,:,1)-latitude_check))==abs(SOCATv2020_grid.latitude(1,:,1)-latitude_check));
lonidx = find(min(abs(SOCATv2020_grid.longitude(:,1,1)-(longitude_check+360)))==abs(SOCATv2020_grid.longitude(:,1,1)-(longitude_check+360)));
latidx=latidx(1);
lonidx=lonidx(1);
% Lat-lon index (2)
latitude_check2 = 45;
longitude_check2 = -130;
latidx2 = find(min(abs(SOCATv2020_grid.latitude(1,:,1)-latitude_check2))==abs(SOCATv2020_grid.latitude(1,:,1)-latitude_check2));
lonidx2 = find(min(abs(SOCATv2020_grid.longitude(:,1,1)-(longitude_check2+360)))==abs(SOCATv2020_grid.longitude(:,1,1)-(longitude_check2+360)));
latidx2=latidx2(1);
lonidx2=lonidx2(1);

figure; hold on;
set(gca,'fontsize',16);
p1=plot(1:12,squeeze(SOCATv2020_grid.pco2_RF_clim(lonidx,latidx,:)),'k-','linewidth',2);
p2=plot(1:12,squeeze(SOCATv2020_grid.pCO2_TA(lonidx,latidx,:)),'r--','linewidth',2);
p3=plot(1:12,squeeze(SOCATv2020_grid.pCO2_DIC(lonidx,latidx,:)),'r:','linewidth',2);
p4=plot(1:12,squeeze(SOCATv2020_grid.pCO2_TMP(lonidx,latidx,:)),'b--','linewidth',2);
p5=plot(1:12,squeeze(SOCATv2020_grid.pCO2_SAL(lonidx,latidx,:)),'b:','linewidth',2);
xlabel('Month of Year');
ylabel('pCO_{2}');
legend({'pCO_{2}' 'pCO_{2(TA)}' 'pCO_{2(DIC)}' 'pCO_{2(Temp.)}' 'pCO_{2(Sal.)}'},...
    'location','northwest','fontsize',16);

figure; hold on;
set(gca,'fontsize',16);
p1=plot(1:12,squeeze(SOCATv2020_grid.pH_clim(lonidx,latidx,:)),'k-','linewidth',2);
p2=plot(1:12,squeeze(SOCATv2020_grid.pH_TA(lonidx,latidx,:)),'r--','linewidth',2);
p3=plot(1:12,squeeze(SOCATv2020_grid.pH_DIC(lonidx,latidx,:)),'r:','linewidth',2);
p4=plot(1:12,squeeze(SOCATv2020_grid.pH_TMP(lonidx,latidx,:)),'b--','linewidth',2);
p5=plot(1:12,squeeze(SOCATv2020_grid.pH_SAL(lonidx,latidx,:)),'b:','linewidth',2);
xlabel('Month of Year');
ylabel('pH');
legend({'pH' 'pH_{(TA)}' 'pH_{(DIC)}' 'pH_{(Temp.)}' 'pH_{(Sal.)}'},...
    'location','northwest','fontsize',16);

figure; hold on;
set(gca,'fontsize',16);
p1=plot(1:12,squeeze(SOCATv2020_grid.OmA_clim(lonidx,latidx,:)),'k-','linewidth',2);
p2=plot(1:12,squeeze(SOCATv2020_grid.OmA_TA(lonidx,latidx,:)),'r--','linewidth',2);
p3=plot(1:12,squeeze(SOCATv2020_grid.OmA_DIC(lonidx,latidx,:)),'r:','linewidth',2);
p4=plot(1:12,squeeze(SOCATv2020_grid.OmA_TMP(lonidx,latidx,:)),'b--','linewidth',2);
p5=plot(1:12,squeeze(SOCATv2020_grid.OmA_SAL(lonidx,latidx,:)),'b:','linewidth',2);
xlabel('Month of Year');
ylabel('\Omega_{A}');
legend({'\Omega_{A}' '\Omega_{A(TA)}' '\Omega_{A(DIC)}' '\Omega_{A(Temp.)}' '\Omega_{A(Sal.)}'},...
    'location','northwest','fontsize',16);

figure; hold on;
set(gca,'fontsize',16);
p1=plot(1:12,squeeze(SOCATv2020_grid.OmA_clim(lonidx2,latidx2,:)),'k-','linewidth',2);
p2=plot(1:12,squeeze(SOCATv2020_grid.OmA_TA(lonidx2,latidx2,:)),'r--','linewidth',2);
p3=plot(1:12,squeeze(SOCATv2020_grid.OmA_DIC(lonidx2,latidx2,:)),'r:','linewidth',2);
p4=plot(1:12,squeeze(SOCATv2020_grid.OmA_TMP(lonidx2,latidx2,:)),'b--','linewidth',2);
p5=plot(1:12,squeeze(SOCATv2020_grid.OmA_SAL(lonidx2,latidx2,:)),'b:','linewidth',2);
xlabel('Month of Year');
ylabel('\Omega_{A}');
legend({'\Omega_{A}' '\Omega_{A(TA)}' '\Omega_{A(DIC)}' '\Omega_{A(Temp.)}' '\Omega_{A(Sal.)}'},...
    'location','northwest','fontsize',16);

%% Latitudinal index
latitude_check = 25;
latidx = SOCATv2020_grid.latitude(1,:,1) < latitude_check;

figure; hold on;
set(gca,'fontsize',16);
p1=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_clim(:,latidx,:),1,'omitnan'),2,'omitnan')),'k-','linewidth',2);
p2=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_TA(:,latidx,:),1,'omitnan'),2,'omitnan')),'r--','linewidth',2);
p3=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_DIC(:,latidx,:),1,'omitnan'),2,'omitnan')),'r:','linewidth',2);
p4=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_TMP(:,latidx,:),1,'omitnan'),2,'omitnan')),'b--','linewidth',2);
p5=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_SAL(:,latidx,:),1,'omitnan'),2,'omitnan')),'b:','linewidth',2);
xlim([1 12]);
xlabel('Month of Year');
ylabel('pH');
ylim([7.98 8.18]);
legend({'pH' 'pH_{(TA)}' 'pH_{(DIC)}' 'pH_{(Temp.)}' 'pH_{(Sal.)}'},...
    'location','northwest','fontsize',16);

figure; hold on;
set(gca,'fontsize',16);
p1=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_clim(:,latidx,:),1,'omitnan'),2,'omitnan')),'k-','linewidth',2);
p2=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_TA(:,latidx,:),1,'omitnan'),2,'omitnan')),'r--','linewidth',2);
p3=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_DIC(:,latidx,:),1,'omitnan'),2,'omitnan')),'r:','linewidth',2);
p4=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_TMP(:,latidx,:),1,'omitnan'),2,'omitnan')),'b--','linewidth',2);
p5=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_SAL(:,latidx,:),1,'omitnan'),2,'omitnan')),'b:','linewidth',2);
xlabel('Month of Year');
ylabel('\Omega_{A}');
xlim([1 12]);
ylim([3.0 3.8]);
legend({'\Omega_{A}' '\Omega_{A(TA)}' '\Omega_{A(DIC)}' '\Omega_{A(Temp.)}' '\Omega_{A(Sal.)}'},...
    'location','northwest','fontsize',16);

%% Distance from shore index
idx = SOCATv2020_grid.distance_from_shore(:,:,1) <= 50 & ...
      ~island(SOCATv2020_grid.lat,SOCATv2020_grid.lon) & ...
      SOCATv2020_grid.latitude(1,:,1) < 50 & ...
      SOCATv2020_grid.latitude(1,:,1) > 22;
idx = repmat(idx,1,1,12);

SOCATv2020_grid.pH_clim(~idx) = NaN;
SOCATv2020_grid.pH_TA(~idx) = NaN;
SOCATv2020_grid.pH_DIC(~idx) = NaN;
SOCATv2020_grid.pH_TMP(~idx) = NaN;
SOCATv2020_grid.pH_SAL(~idx) = NaN;
figure; hold on;
set(gca,'fontsize',16);
p1=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_clim,1,'omitnan'),2,'omitnan')),'k-','linewidth',2);
p2=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_TA,1,'omitnan'),2,'omitnan')),'r--','linewidth',2);
p3=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_DIC,1,'omitnan'),2,'omitnan')),'r:','linewidth',2);
p4=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_TMP,1,'omitnan'),2,'omitnan')),'b--','linewidth',2);
p5=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.pH_SAL,1,'omitnan'),2,'omitnan')),'b:','linewidth',2);
xlim([1 12]);
xlabel('Month of Year');
ylabel('pH');
ylim([7.95 8.15]);
legend({'pH' 'pH_{(TA)}' 'pH_{(DIC)}' 'pH_{(Temp.)}' 'pH_{(Sal.)}'},...
    'location','northwest','fontsize',16);

SOCATv2020_grid.OmA_clim(~idx) = NaN;
SOCATv2020_grid.OmA_TA(~idx) = NaN;
SOCATv2020_grid.OmA_DIC(~idx) = NaN;
SOCATv2020_grid.OmA_TMP(~idx) = NaN;
SOCATv2020_grid.OmA_SAL(~idx) = NaN;
figure; hold on;
set(gca,'fontsize',16);
p1=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_clim,1,'omitnan'),2,'omitnan')),'k-','linewidth',2);
p2=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_TA,1,'omitnan'),2,'omitnan')),'r--','linewidth',2);
p3=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_DIC,1,'omitnan'),2,'omitnan')),'r:','linewidth',2);
p4=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_TMP,1,'omitnan'),2,'omitnan')),'b--','linewidth',2);
p5=plot(1:12,squeeze(mean(mean(SOCATv2020_grid.OmA_SAL,1,'omitnan'),2,'omitnan')),'b:','linewidth',2);
xlabel('Month of Year');
ylabel('\Omega_{A}');
xlim([1 12]);
ylim([2.1 2.9]);
legend({'\Omega_{A}' '\Omega_{A(TA)}' '\Omega_{A(DIC)}' '\Omega_{A(Temp.)}' '\Omega_{A(Sal.)}'},...
    'location','northwest','fontsize',16);