function [mData] = d2ms(Data,month_length)
%This function converts the daily time series 'Data' to monthly time series
%by accumulating the values
% mData is the resulting monthly time series
% Data is the original daily data
% month_length is the array of the number of days of each month (e.g. for the period 1981-1982 N = [31 28 31 30 31
% 30 31 31 30 31 30 31 31 28 31 30 31
% 30 31 31 30 31 30 31] 

% By Mosisa Tujuba Wakjira
% Version 1: 03 August 2022

cumML = cumsum(month_length);
mData = zeros(length(month_length),1);
for i=1:length(month_length)
    mData(i)=sum(Data(cumML(i)-(month_length(i)-1):cumML(i)),"omitmissing");
end
end