
%% Saxton and Rawls (2006) PTF to compute soil moisture constants 
% Permanent wilting point (PWP)
% Author: Mosisa Tujuba Wakjira
 wp0 = -(0.024*fsand)+(0.487*fclay)+(0.006*fsom+(0.005*(fsand.*fsom)))-...
     (0.013*(fclay.*fsom))+(0.068*(fsand.*fclay))+0.031;
fWP = wp0+(0.14*wp0)-0.02;
 
 % Field capacity (FC) and saturation water content 
 
 fc0 = (-0.251*fsand)+(0.195*fclay)+(0.011*fsom)+(0.006*(fsand.*...
        fsom))-(0.027*(fclay.*fsom))+(0.452*(fsand.*fclay))+0.299;
    
 fFC = fc0 + ((1.283*(fc0.^2))-(0.374*fc0)-0.015);
 
 st0 = (0.278*fsand)+(0.034*fclay)+(0.022*fsom)-(0.018*(fsand.*...
        fsom))-(0.027*(fclay.*fsom))-(0.584*(fsand.*fclay))+0.078;
    
 st1 = st0+((0.636*st0)-0.107);
    
 fSat = (fFC+st1)+((-0.097*fsand)+0.043);
 clear wp0 fc0 st0 st1