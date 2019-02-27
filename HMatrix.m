function C = HMatrix(H)
% ***************************************************************************
% Description: return a matrix of h 
%
% Parameters:
% h: an array of unit holding cost in each period
%
%
% author: Zhen Chen
% time: 2019-02-25, 17:16
% ***************************************************************************

d = length(H);
tempH = repmat(H', 1, d);
tempH = triu(tempH, 1); 
C = cumsum(tempH); 
end