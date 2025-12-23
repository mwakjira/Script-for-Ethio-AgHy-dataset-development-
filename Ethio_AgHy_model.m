%%%%%%%%%%%% SOIL WATER BALANCE SIMULATION MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%% FOR THE ETHIOPIAN RFA REGION %%%%%%%%%%%%%%%%%
% Author: Mosisa Tujuba Wakjira

% This water balance model simulates hydrological fluxes at the surface and within the top
% 60cm of the root zone. The rainfall-runoff partitioning is performed using curve number method.

%% Import data
clc
clear
cd('D:\AgHyCA\Waterlogging');
load('WLR_500m_input.mat')


%% Parameters
% --- Root zone depth and moisture depletion factor ---
Zr = 0.60; % root zone depth of the model (single layer)
f = 0.55; % available moisture depletion factor, which depends on crop type 
% and ETc (Table 22, Allen et al., 1998). The value assumed here
% corresponds to four major cereal crops in Ethiopia (maize, teff, sorghum and wheat)

% ---Curve number ---
% The CN used in this study corresponds to agricultural land use, and
% specific to the hydrologic soil groups. The fixed curve number (fCN) was adjusted for slope
% according Williams et al., 2023 (see script 'WLR_data_preprocessing.m').
% The average slope-adjusted fCN was also computed for dry and wet
% conditions and used as reference values to calculate the dynamic (daily)
% CN corresponding to the daily moisture conditions

% CN for wet and dry AMC
CNw = fCN.*exp(0.00673*(100-fCN)); % CN for wet soils (above FC)
CNd = fCN-(20*(100-fCN))./((100-fCN)+exp(2.533-0.0636*(100-fCN))); % for dry soils (below WP)
fFW = 0.5*(fFC+fWP); % intermediate point between FC and WP (fFW) corresponding the average CN

% --- Soil water retention parameters (theta) ---
thSat = 1000*fSat*Zr; % soil moisture at saturation in mm 
thFC = 1000*fFC*Zr; % soil moisture at field capacity in mm
thWP = 1000*fWP*Zr; % soil moisture at wilting point in mm
thRes = 0.01*(0.145 + 0.253*(100*0.58*fsom) + 0.033*(100*fclay)); % residual moisture content (based on Poeplau et al., 2015)
fRAW = 1000*f*(fFC-fWP)*Zr; % fixed readily available water (RAW)

% --- variables to be stored during the simulation ---
Ro = nan(10957,1); % runoff time series 
Sm = nan(10957,1); % soil moisture 
et = nan(10957,1); % actual ET
Ad = nan(10957,1); % aeration deficit
% d = nan(10957,1); % deep percolation
RO = nan(2120,2250,360); 
SMv = nan(2120,2250,360);
ETa = nan(2120,2250,360);
% Dp = nan(2120,2250,360);
NoE = nan(2120,2250,30); % Number of waterlogging events
ROD = nan(2120,2250,30); % average magnitude (Relative Oxygen Deficit)
DoE = nan(2120,2250,30); % total duration of waterlogging 
% EOD = nan(2120,2250,360); % Earliest onset date of waterlogging
% LOD = nan(2120,2250,360); % Latest onset of date of waterlogging

% c = 0;

% --- Simulations at each grid ---
for i=0:211 % Rainfall and ETo are available only at 5km resolution (10x coarser)
    for j=0:224
        for m=1:10
            for n=1:10
            ii = 10*i+m; jj = 10*j+n;
            if isnan(rfa_HR(ii,jj))
                continue
            end
        
            %% Initial condition and model stablization period
            % As an intial condition, we assume that soil moisture on Jan
            % 1, 1981 is 20% higher than moisture content at WP. We
            % provided a 3-year stablization period for the model, i.e.,
            % the model was run for three years with the initial condition
            % provided. 

            r = squeeze(RF(i+1, j+1, :)); P = [r(1:1095); r]; % daily RF
            p = squeeze(pet(i+1, j+1, :)); ET = [p(1:1095); p]; % daily ETo
            
            sm0 = 1.2*fWP(ii,jj); % We set an arbitrary initial condition 
           
            for k=1:12052 % out of which 1138 days are stablization period
                sm = 1000*Zr*(sm0-thRes(ii,jj)); % soil moisture in mm. Soil 
                % moisture was allowed to deplete below the wilting point 
                % down to residual moisture, accounting for soil 
                % evaporation and vapor losses (avoids static "dead
                % volume"). This assumption follows the fact that soil drying (evaporation)
                % continues (Seneviratne et al., 2010) beyond WP although this is 
                % process has no biological significance plant uptake. It also follows the 
                % widely used soil moisture parameterization of soil water content 
                % by van Genuchten (1980), which allows soil water content to smoothly 
                % toward residual moisture content without cutoff at WP.
        
                Storage = thSat(ii,jj) - sm; %available storage (mm)

                % --- Adjust CN for antecedent moisture conditions ---
                if sm0 >= fFC(ii,jj)
                    CN = CNw(ii,jj);
                elseif sm0 <= fWP(ii,jj)
                    CN = CNd(ii,jj);
                elseif sm0 > fWP(ii,jj) && sm0 <= fFW(ii,jj)
                    CN = CNd(ii,jj)*(fFW(ii,jj)-sm0)/(fFW(ii,jj)-...
                        fWP(ii,jj))+fCN(ii,jj)*(1-(fFW(ii,jj)-sm0)/...
                        (fFW(ii,jj)-fWP(ii,jj)));
                else
                    CN = fCN(ii,jj)*(sm0-fFW(ii,jj))/(fFC(ii,jj)-...
                        fFW(ii,jj))+CNw(ii,jj)*(1-(sm0-fFW(ii,jj))/...
                        (fFC(ii,jj)-fFW(ii,jj)));
                end

                % --- calculate maximum surface retention S ---
                S = 254*(100/(0.96*CN)-1); % 0.96 is a pre-determined 
                % scaling factor optimized using observed runoff 
                
                % --- calculate surface runoff Q ---
                if P(k)>0.071*S
                    Q = (P(k) - 0.071*S)^2/(P(k)+1.349*S);
                    Ip = P(k) - Q; % possible infiltration
                else
                    Q = 0;
                    Ip = P(k);
                end

                % We assume that runoff on flat slope of less than 2%
                % (FAO, 2006 guidelines for soil description), not suffiently drain and hence will have
                % high infiltration opportunity time due to the possibility
                % of ponding. Thus, we assumed that runoff is a linear
                % function slope such that Q=0 at slope=0
                if P(k)>0.071*S && slope(ii,jj) <= 2
                   Q = (slope(ii,jj)/2)*((P(k) - 0.071*S)^2/(P(k)+1.349*S));
                   Ip = P(k) - Q; 
                end
                % Maximum storage capacity of soil is equal to moisture 
                % content at saturation (see Raes et al., 2012, chapter 3).
                % Therefore the dynamic soil stroage is ranges from 0
                % (saturated) to 1000*(fSat - thRes)*Zr 
                if Ip>Storage 
                   Q = Q + (Ip - Storage);
                   sm = sm + Storage; % antecedent sm and infiltration
                else
                   sm = sm + Ip;
                end

                % --- calculate water stress coefficient Ks [-] and ETa ---
                theta = sm/(1000*Zr); % sm in mm
                if theta < 0.5*fWP(ii,jj) 
                    Ks = 0;
                    % The use of 0.5*WP follows the above assumption regarding soil drying.
                    % Transpiration stops at WP but drying continues (Seneviratne et al., 2010). 
                    % We partly represent the water-filled pore dynamics due to extreme drying 
                    % by assuming 0.5*WP as the lowest margin. This reflects extreme dry-end 
                    % conditions and users may adjust or introduce a parameter to calibrate it
                    % depending on application. 
                elseif sm <= f*(thFC(ii,jj) - thWP(ii,jj))
                    Ks = sm/(f*(thFC(ii,jj) - thWP(ii,jj)));
                else
                    Ks = 1;
                end
                ea = Ks * ET(k);
                sm = sm - ea;
            
                % --- Calculate deep percolation, Dp [mm] ---
                % Dp is limited by soil hydraulic properties, which is
                % mainly defined by soil textural composition and organic
                % matter content. 
                theta = sm/(1000*Zr);
                sfc = fFC(ii,jj); scl = fclay(ii,jj); ssn = fsand(ii,jj);
                sst = thRes(ii,jj); srw = fRAW(ii,jj);
                
                        if theta > sfc && scl > 0.40 &&...
                           ssn < 0.45 && 1-(scl + ssn)<0.4 % if the soil is saturated and clay
                           Kt = 10*60*24*(0.231*(theta-sst)^7.6512); % unsaturated hydraulic conductivity (based on Zayani et al., 1992), converted from cm/min to mm/day
                           dp = Kt;
                        elseif theta > sfc % if the soil is saturated and clay
                           dp = sm - thFC(ii,jj);   
                        else
                           dp = 0;
                        end
                sm = sm - dp;
                sm0 = sm/(1000*Zr)+thRes(ii,jj); % update the soil 
                % volumetric soil moisture to serve as antecedent soil 
                % moisture for the next day
            
                %% store the relevant water balance components
                if k>1095
                   Ro(k-1095) = Q*sm/sm;
                   Sm(k-1095) = sm0;
                   et(k-1095) = ea;
                   % d(k-1095) = dp;
                   s = NaN;
                   if sm0 > sfc
                      s = (sm0-sfc)/(fSat(ii,jj)-sfc);  
                   end
                   Ad(k-1095) = s; % soil aeration deficit
                end
            end

            % --- Waterlogging conditions
            leapY = [0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0];
            leapY = leapY == 1;
            for kk = 0:29 % during the period 1981-2010
                if ~leapY(kk+1)
                    sm = Sm(365*kk+1:365*kk+365);
                else
                    sm = Sm(365*kk+1:365*kk+366);
                end
                sme = sm-fFC(ii,jj); % excess soil moisture in %
                % indentify the days of excess water
                indx = double(sme>0); % index for the days of the year with excess soil moisture
                if isempty(indx) 
                   continue
                end
                for s = 1:length(indx)-2
                    if indx(s)~=indx(s+1) && indx(s)==indx(s+2) 
                    % a non-waterlogged day between multiple waterlogged day was considered also a waterlogged day, or vice versa.
                    indx(s+1) = indx(s);
                    end
                end

                indx = find(indx==1); % days with aeration deficit
                indn1 = find(diff(indx)>1); 

               
                if ~isempty(indx)
                   NoE(ii,jj,kk+1) = length(indn1)+1; % Number of waterlogging events 
                   DoE(ii,jj,kk+1) = length(indx); % total duration of waterlogging per year
                   ROD(ii,jj,kk+1) = 100*mean(sme(indx))/(fSat(ii,jj)-fFC(ii,jj)); % average magnitude (Relative Oxygen Deficit as a percent air filled porosity of the soil)
                   % EOD(ii,jj,kk+1) = indx(1); % Earliest onset date of waterlogging
                   % LOD(ii,jj,kk+1) = indx(end); % Latest onset of date of waterlogging
                end
                clear indn1 indn2
            end

             %% Runoff and erosion metrics
             % Monthly soil moisture and runoff
             RO(ii,jj,:) = d2ms(Ro,monthL);
             SMv(ii,jj,:) = d2mm(Sm,monthL);
             ETa(ii,jj,:) = d2ms(et,monthL);
             % Dp(ii,jj,:) = d2ms(dp,monthL);
            end
        end
    end
    % c = c+1
end

% --- calculate waterlogging intensity ---
WLI = nan(2120,2250,30);
for i = 1:30
    WLI(:,:,i) = 10*0.6*ROD(:,:,i).*(fSat-fFC).*DoE(:,:,i)/600; % ROD converted to fraction, hence time 10, not 1000.
end

% --- save the output agrohydrological variables
save('WLR_500m_output.mat','RO','ETa','SMv','ROD','NoE','DoE','WLI','-v7.3');
clear cumML clear p m n sm sm0 ii jj P ET S CN kk Q Q_adj Ip I Storage theta theta_r c i j sm_v ea dp Kt Ks y yy xx scl ssn sst ssl sfc r





