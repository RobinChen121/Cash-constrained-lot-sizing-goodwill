function [finalXMn, finalB, satDemand, finalIniB, finalIniW, finalW] = ComputeMkn(xMn, lastX, thisIniB, ... 
                                    thisIniW, lastB, lastW, lastSatDemand, d, p, c, s, h, beta, needLossless, thisLoanLength, ratePay)
                                
% ***************************************************************************
% Description: get objective matrix, constraint matrix for an ordering
% round to compute the maximum cash increment
%
% Parameters:
% xMn: ordering plan in an ordering cycle from period m to period n
% lastX: an array, previous optimal ordering plan
% thisIniB: an array, initial cash for each ordering cycle in this ordering round
% thisIniW: an array, initial lost sale for each ordering cycle in this ordering round
% lastB: end-of period cash in the final period of previous ordering plan
% lastW: lost sale quantity in the final period of previous ordering plan
% lastSatDemand: an array, satifying quantity of the previous ordering plan
% d: (1*n) demands
% p: (1*n) prices
% s: (1*n) fixed ordering costs
% c: (1*n) unit vari ordering costs
% h: (1*n) unit holding costs
% beta: goodwill loss rate
% needLossless: binary variable, whether need to add the less lost sale
% constraint
% thisLoanLength: loaning length in this ordering round 
% ratePay: total loan and interest needed to pay back
%
% Decision variables:
% finalXMn: final ordering plan
% finalB: final end-of-period cash
% satDemand: final satisfying quantity
% finalIniB: final initial cash
% finalIniW: final initial loss sale
% finalW: final lost sale in the final period of the ordering round
%
% author: Zhen Chen
% time: 2019-02-19, 17:16
% ***************************************************************************                                
                                
%% initial setting, the ordering round: from period m to period n, include two ordering cycles: m ~ k, k + 1 ~ n
finalB = lastB;
options = optimset('TolFun', 1e-6, 'MaxIter', 20);
n = length(p);
m = 1;
iniB = thisIniB(1);
iniW = thisIniW(1);
finalXMn = xMn;
[~, orderIndex] = find(xMn == 1);                                
                                
%% calling ComputefAb to get f A b
[f, A, b] = ComputefAb(xMn, iniB, iniW, lastW, d(m : n), p(m : n), c(m : n), s(m : n), h(m : n), beta, zeros(1, n - m + 1),...
                       needLossless, thisLoanLength, ratePay);
lb = zeros(1,n); ub = d(1 : n);

%% calling linprog to compute satDemand
[satDemand, fval, exitflag] = linprog(f, A, b, [], [], lb, ub, [], options);

% if no solution, relax the effective demand constraint
if exitflag ~= 1
    [satDemand, ~, exitflag1] = linprog(f, A(1 : n - 1 + sum(xMn), :), b(1 : n - 1 + sum(xMn), :), [], [], lb, ub, [], options);  
    if exitflag1 ~= 1
        satDemand = zeros(1, n - m + 1)';
    end
end

%% check whether cash is increasing or not, and update
orderLastLength = diff(orderIndex);
orderLastLength(end + 1) = n - m - orderIndex(end) + 2; 
finalIniW = zeros(1, length(thisIniW));
finalIniB = zeros(1, length(thisIniW));
if needLossless == 0 && thisLoanLength == n % need to pay back loan in period n
    lastB = lastB - ratePay; 
end
if exitflag ~= 1 || finalB < lastB - 1e-2
    finalB = lastB;
    finalIniB = thisIniB;
    finalIniW = thisIniW;
    if needLossless == 0 % two situations
        satDemand = [lastSatDemand; 0]; finalXMn = [lastX,0];
    else
        satDemand = lastSatDemand; finalXMn = lastX;
    end 
else
    finalIniB(1) = iniB;
    % recompute finalIniW when solution is like 0, 0, 0, 1
    if orderIndex(1) > 1
        [~,~,wFlow] = ComputeEd(zeros(1, orderIndex(1) - 1), d(1 : orderIndex(1) - 1), iniW, beta);
        finalIniW(1) = wFlow(end);
    else
        finalIniW(1) = iniW; 
    end
    for xNum = 1 : sum(xMn) - 1
        kStart = orderIndex(xNum); % starting period index of cycle m~k
        kLast = kStart + orderLastLength(xNum) - 1;
        finalIniB(xNum + 1) = -A(kLast,:) * satDemand - s(1 : kLast) * xMn(1 : kLast)' + iniB; 
        if thisLoanLength > 0 && orderIndex(xNum + 1) > thisLoanLength % need to pay back loan in this situation
            finalIniB(xNum + 1) = finalIniB(xNum+1) - ratePay;
        end
    end 
end
[~, ~, w] = ComputeEd(satDemand, d, iniW, beta);
finalW = w(n);
for xNum = 1 : sum(xMn) - 1
    kStart = orderIndex(xNum);
    kLast = kStart + orderLastLength(xNum) - 1;
    finalIniW(xNum + 1) = w(kLast);
end

end