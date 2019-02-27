function [f, A, b] = ComputefAb(x, iniB, iniW, preFinalW, d, p, c, s, h, beta, isLossMuch, needLossless, thisLoanLength, ratePay) 

% ***************************************************************************
% Description: get objective matrix, constraint matrix for an ordering
% round to compute the maximum cash increment
%
% Parameters:
% x: the ordering plan in the ordering round
% iniB: initial cash in the ordering round
% iniW: initial lost sale in the ordering round
% preFinalW: final lost sale in the last peirod of previous plan
% WW: end-of period lost sale in each ordering cycle of the matrix
% n: length of x
% d: (1*n) demands
% p: (1*n) prices
% s: (1*n) fixed ordering costs
% c: (1*n) unit vari ordering costs
% h: (1*n) unit holding costs
% beta: goodwill loss rate
% isLossMuch: (1*n) binary variables, whether this period lost sale too
% much
% needLossless: binary variable, whether need to add the less lost sale
% constraint
% thisLoanLength: loaning length in this ordering round 
% ratePay: total loan and interest needed to pay back
%
% Decision variables:
% f: coefficient vector for objective function
% A: constraint matrix left side
% b: constraint matrix right side
% author: Zhen Chen
% time: 2019-02-19, 17:16
% ***************************************************************************

n = length(p);  % length from period i to period j
xSum = (sum(x)); % total ordering launching times
[~, orderingIndex] = find(x == 1);  % record the ordering index 
orderLastLength = diff(orderingIndex); % record ordering cycle lengthes in the ordering round
orderLastLength(end + 1) = n - orderingIndex(end) + 1;

%% compute f
tempH = zeros(n, n); tempC = zeros(n, n); tempS = zeros(1, n);
for i = 1 : xSum
    m = orderingIndex(i);
    t = orderingIndex(i) + orderLastLength(i) - 1; % period m to period t is an ordering cycle
    tempH(m : t, m : t) = HMatrix(h(m : t));
    tempH(t : n, :) = repmat(tempH(t, :), n - t + 1, 1); tempC(m : n, m : t) = c(m);
    if i == 1
        tempS(m : t) = s(1);
    else
        tempS(m : t) = s(m) + tempS(m - 1);
    end
end
f = tempC(n, :) - p(1 : n) + tempH(n, :); 
f = f';

%% compute A
%% the first n - 1 constraint is cash flow constraint, end-of-period cash
% should be higher than 0
if needLossless == 1
    A = zeros(n + xSum + n, n);
    b = zeros(n + xSum + n, 1); % add one more constraint 
else
    A = zeros(n + xSum + n - 1, n);
    b = zeros(n + xSum + n - 1, 1);
end
tempP = tril(repmat(p(1 : n), n, 1)); 
A(1 : n - 1, :) = tempC(1 : n - 1, :) - tempP(1 : n - 1, :) + tempH(1 : n - 1, :);
b(1 : n - 1) = iniB - tempS(1 : n - 1);
if thisLoanLength > 0
    b(thisLoanLength : n - 1) = b(thisLoanLength : n - 1) - ratePay;
end

%% ordering quantity constrint, from the n th to the n + xSum - 1 th
for i = 1 : xSum
    m = orderingIndex(i); mLastLength = orderLastLength(i);
    if i==1
        A(n, m : m + mLastLength - 1) = c(1);
        b(n) = iniB - s(1);
    else
        A(n + i - 1, :) = A(m - 1, :) + [zeros(1, m - 1), c(m) * ones(1, mLastLength), zeros(1, n - m - mLastLength + 1)];
        b(n + i - 1, :) = b(m - 1) - s(m);
    end
end

%% effective demand constraint: v_t <= Ed_t, from n + xSum th to n + xSum + n - 1 th
Ed = d;
if orderingIndex(1)>1
    [Ed(1 : orderingIndex(1) - 1), ~, wFlow] = ComputeEd(zeros(1, orderingIndex(1) - 1), d(1 : orderingIndex(1) - 1), iniW, beta);
    iniW = wFlow(end);
end
for i = 1 : orderingIndex(1) - 1
    A(n + xSum + i - 1, i) = 1; b(n + xSum + i - 1, i) = 0; % v_t = 0 at this situation
end

Ed(orderingIndex(1)) = max(0, d(orderingIndex(1)) - beta * iniW);
preIndexLossMuch = orderingIndex(1) - 1;
i = orderingIndex(1);
while i <= n
    if i + 1 <= n 
        rowIndex = n + xSum + i -1;
        if isLossMuch(i + 1) ~= 1
            tempBetaArray = BetaArray(i - preIndexLossMuch, beta);
            A(rowIndex, preIndexLossMuch + 1 : i) = tempBetaArray; 
            b(rowIndex) = Ed(preIndexLossMuch + 1 : i) * tempBetaArray';
        else
            tempBetaArray = BeitaZhengfu(i - preIndexLossMuch, beita);
            A(rowIndex, preIndexLossMuch + 1 : i) = tempBetaArray;
            b(rowIndex) = Ed(preIndexLossMuch + 1 : i) * tempBetaArray' - 2 * d(i + 1); % effective demand constraint for period i
            A(rowIndex + 1, i + 1) = 1; % effective demand constraint for period i + 1
            preIndexLossMuch = i + 1;
            i = i + 1;
        end
    end
    i=i+1;
end
if isLossMuch(n) ~= 1
    tempBetaArray = BetaArray(n - preIndexLossMuch, beta);
    A(n + xSum + n - 1, preIndexLossMuch + 1 : n) = tempBetaArray; b(n + xSum + n - 1) = Ed(preIndexLossMuch + 1 : n) * tempBetaArray';
end

%% one constraint: can not lose more than the previous plan
if needLossless == 1&& n - preIndexLossMuch > 0
    tempBetaArray = BetaArray(n - preIndexLossMuch, beta);
    A(n + xSum + n, preIndexLossMuch+1 : n) = -tempBetaArray;
    b(n + xSum + n) = preFinalW - Ed(preIndexLossMuch + 1: n) * tempBetaArray'+1e-2;
end

end