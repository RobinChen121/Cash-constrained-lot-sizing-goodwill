function [Ed, isLossTooMuch, w] = ComputeEd(satDemand, d, iniW, beta)

% ***************************************************************************
% Description: compute effective demand
%
% Parameters:
% N: length of satDemand
% satDemand: (1*N) satisfying demand in each period
% d: (1*N) demands
% iniW: lost sale quantity at the beginning 
% beta: goodwill loss rate
%
%
% Decision variables:
% Ed: (1*N) effective demand in each period
% isLossTooMuch: (1*N) binary variable judging whether this period loss too
% much
% w: (1*N) demand shortage (lost sales) in each period
%
% author: Zhen Chen
% time: 2019-02-19, 22:16
% ***************************************************************************


N = length(d); 
Ed = zeros(1, N);
isLossTooMuch = zeros(1, N);
w = Ed;
Ed(1) = max(0, d(1) - beta * iniW);
for i = 1 : N
    if i > 1
        if d(i) - beta * w(i - 1) > 1e-2
            Ed(i) = d(i) -beta * w(i - 1);
        else
            isLossTooMuch(i) = 1; 
            Ed(i) = 0;
        end
    end
    w(i) = Ed(i) - satDemand(i);
end

end