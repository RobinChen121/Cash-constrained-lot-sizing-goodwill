function [finalIniB, finalIniW] = UpdateIniBW(xMn, IniB, IniW, iniBMn, iniWMn)

% ***************************************************************************
% Description: update the initial cash B and initial lost sale quantity w in
% the each ordering cycle in an ordering round
%
% Parameters:
% xMn: an ordering round from period m to period n
% iniB: an array, initial cash of each ordering cycle
% iniW: an array, initial lost sale of each ordering cycle
% BB: a matrix, end-of period cash in each ordering cycle of the matrix
% WW: a matrix, end-of period lost sale in each ordering cycle of the matrix
%
%
% Decision variables:
% finalIniB: an array, initial cash for each ordering cycle in the ordering round
% finalIniW: an array, initial lost sale for each ordering cycle in the ordering round
%
% author: Zhen Chen
% time: 2019-02-19, 16:16
% ***************************************************************************

finalIniB = iniBMn; finalIniW = iniWMn;
[~, orderingIndex] = find(xMn == 1);
n = length(xMn);
if ~isempty(orderingIndex)
    orderingLastLength = diff(orderingIndex);
    orderingLastLength(end + 1) = n - orderingIndex(end) + 1; 
end
for xNum = 1 : sum(xMn)
    t1 = orderingIndex(xNum);
    t2 = orderingIndex(xNum) + orderingLastLength(xNum) - 1;
    finalIniB(t1, t2) = IniB(xNum);
    finalIniW(t1, t2) = IniW(xNum);
end

end