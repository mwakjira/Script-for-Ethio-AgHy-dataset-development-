%%%%%%% Script to process topographic, soil and Curve Number data %%%%%%%
% The processed data were used as input for the agrohydrological modeling
% framework for the Ethio-AgHy dataset generation and waterloggging impact
% modeling across the RFA region of Ethiopia
% Last update: 5 November 2025
% By Mosisa Tujuba Wakjira, mosisatujuba@gmail.com


%% --- PROCESS THE SLOPE DATA ---
% We use the shuttle radar topographic mission (SRTM) 30 DEM data to
% copmpute slope. Data last accessed on 21.03.2025

cd('D:\AgHyCA\Waterlogging\rasters_SRTMGL1');
d = dir;

% Get the DEM data
dem30 = [];
for i = 3:6
    fname = d(i).name;
    xx = double(imread(fname,'tif'));
    if ismember(i,[3,4,5])
        xx(:,end) = [];
    end
    dem30 = [dem30,xx];
end
dem30(dem30<0) = NaN; 
clear x xx

% Calculate 30m slope (%)
dx = 30; dy = 30;
[dz_dx, dz_dy] = gradient(dem30, dx, dy);
slop30 = sqrt(dz_dx.^2 + dz_dy.^2) * 100;
slop30(slop30>500) = 500;
clear dx dy dz_dy dz_dx

% Calculate 500m slope (%)
nlt=linspace(15.0463,4.45559,2120);
nln=linspace(32.8525,44.0855,2250);
[nln,nlt]=meshgrid(nln,nlt);
lat = linspace(15.0463,4.45559,size(slop30,1))';
lon = linspace(32.8525,44.0855,size(slop30,2))';
slope = interp2(lat,lon,slop30',nlt,nln,'nearest');

cd('D:\AgHyCA\Waterlogging')
save('wlr_topodata.mat','slop30','dem30','slope','-v7.3');


%% --- PROCESS THE CURVE NUMBER DATA ---
% The curve number (CN) data used in this analysis is derived from the USDA
% lookup table (USDA,1985) considering the hydrologic soil group dataset
% (HYSOGs250m, Ross et al., 2018, accessed on 28.03.2025).

% The the gridded CN considers croplands dominated by small grains and legumes
% mainly broadcasted in what is approximately contour farming. The average
% of good and poor conditions were considered to account for the variations
% in the agricultural land use conditions. Bare soil conditions have also been considered 
% in order to represent the soil surface conditions at the fallow season, and early growing
% stage after sowing.

% Load the raw data
cd('D:\AgHyCA\Waterlogging');
CNA2=nan(4240,4500); % average CN for agricultural land (CNA2) from the USDA HEH lookup table (USDA, 1985)
HSG = double(imread("HSG250.tif","tif")); % soil hydrologic group 
HSG(HSG>10) = NaN;
CNA2(HSG==1)=64;
CNA2(HSG==2)=75;
CNA2(HSG==3)=83;
CNA2(HSG==4)=86;

%% Adjusting CN for slope

% The tabulated CN is primarily suitable for slope conditions of < 5%,
% which entails adjustment of CN for use in a wider range of slope
% conditions. 

% We make the adjustment to slope as proposed by Williams et al., 2023 as
% implemented in APEX model. 

S = 254*(100./CNA2-1); % calculate the maximum surface retention S using CNA2
S_adj = S.*(1.1-slope./(slope+2.718.^(3.7+0.02117*slope)));
fCN = 100./(S_adj/(254)+1); % fixed, slope adjusted CN

% lat and lon for the high resolution agrohydrological modelling to map
% waterlogging and erosion risk
rfa_HR = double(imread("rfa_HR.tif","tif")); rfa_HR(rfa_HR<0)=NaN; % RFA grids
wlr_Lat = linspace(15.0463,4.45559,2120)';
wlr_Lon = linspace(32.8525,44.0855,2250)';

save("WLR_curve_number","HSG", "CNA2","wlr_Lon","wlr_Lon","fCN");

%% This code collects and preprocesses the data required for the Ethio-AgHy and waterlogging mapping projects
%% Soil parameters (250 x 250m)
% Load the raw data
cd('D:\PhD\soil data'); % The soilgrids data collected during the PhD was used
load("soil_parameters.mat","fclay","fsand","fsom","rfa"); % soilgrids data

% make grid alignment
fclay(4785:end,:)=[]; fclay(1:399,:)=[]; fclay(:,5342:end)=[]; fclay(:,1:377)=[];
fsand(4785:end,:)=[]; fsand(1:399,:)=[]; fsand(:,5342:end)=[]; fsand(:,1:377)=[];
fsom(4785:end,:)=[]; fsom(1:399,:)=[]; fsom(:,5342:end)=[]; fsom(:,1:377)=[];

% make grid size alignment
nlt=linspace(15.0463,4.45559,2120);
nln=linspace(32.8525,44.0855,2250);
[nln,nlt]=meshgrid(nln,nlt);
lat = linspace(15.0463,4.45559,size(fclay,1))';
lon = linspace(32.8525,44.0855,size(fclay,2))';
fclay = interp2(lat,lon,fclay',nlt,nln,'linear');
fsand = interp2(lat,lon,fsand',nlt,nln,'linear');
fsom = interp2(lat,lon,fsom',nlt,nln,'linear');

% Compute the soil moisture constants using PTF (Saxton and Rawls, 2003)
cd("D:\AgHyCA\Agrohydrological elaboration\Ethio_AgHy");
run Ethio_AgHy_PTF_SoilMoistureRetentionPars.m


%% --- SOIL TEMPERATURE CLIMATOLOGY DATA WATERLOGGING IMPACT ASSESSMENT ---

fTsoil = (double(imread("soilT0_5.tif","tif")) +...
    double(imread("soilT5_15.tif","tif")))/2;
fTsoil(fTsoil<0) = NaN;

cd("D:\AgHyCA\Agrohydrological elaboration\Ethio_AgHy");
save('WLR_inputParameters_250m.mat','wlr_Lat','wlr_Lon','fCN',...
    'fsom','fsand','fclay','fFC','fWP','fSat','slope','rfa_HR',...
    'monthL','rfa_shp','fTsoil','pet','RF','-v7.3');