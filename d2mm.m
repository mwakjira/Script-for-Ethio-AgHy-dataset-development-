function [mData] = d2mm(Data,month_length)
%This function converts the daily time series 'Data' to monthly time series
%by taking the mean of the values
% mData is the resulting monthly time series
% Data is the original daily data
% month_length is the array of the number of days of each month (e.g. for the period 1981-1982 N = [31 28 31 30 31
% 30 31 31 30 31 30 31 31 28 31 30 31
% 30 31 31 30 31 30 31] 

% Version 1: 03 August 2022
% Author: Mosisa Tujuba Wakjira

cumML = cumsum(month_length);
mData = zeros(length(month_length),1);
for ii=1:length(month_length)
    mData(ii)=mean(Data(cumML(ii)-(month_length(ii)-1):cumML(ii)),"omitmissing");
end
end

