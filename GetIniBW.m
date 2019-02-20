function [iniB,iniW] = GetIniBW(xMn, B0, w0, BB, WW)

% ***************************************************************************
% Description: get the initial cash B and initial lost sale quantity w in
% the each ordering cycle in an ordering round
%
% Parameters:
% xMn: an ordering round from period m to period n
% B0: initial cash of the ordering ournd
% W0: initial lost sale quantity of the ordering round
% BB: end-of period cash in each ordering cycle of the matrix
% WW: end-of period lost sale in each ordering cycle of the matrix
%
%
% Decision variables:
% iniB: an array, initial cash for each ordering cycle in the ordering round
% iniW: an array, initial lost sale for each ordering cycle in the ordering round
%
% author: Zhen Chen
% time: 2019-02-19, 16:16
% ***************************************************************************

[~, orderingIndex] = find(xMn == 1); n = length(xMn);
orderLastLength = diff(orderingIndex); orderLastLength(end + 1) = n - orderingIndex(end) + 1; %记录启动开始阶段，持续时长
iniB = zeros(1, sum(xMn)); iniW = zeros(1, sum(xMn));
iniB(1) = B0; iniW(1) = w0; 
for xNum = 1 : sum(xMn) - 1
    t1 = orderingIndex(xNum); % ordering cycle starting index
    t2 = orderingIndex(xNum) + orderLastLength(xNum) - 1; % ordering cycle ending index
    iniB(xNum+1) = BB(t1, t2);
    iniW(xNum+1) = WW(t1, t2);
end

end