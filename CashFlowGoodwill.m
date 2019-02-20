function [x, y, w, I, B, whetherAdjustPlan, whetherMoveOrderQuantity] = CashFlowGoodwill(d, p, s, c, h, beta, B0, BL, TL, rL)

% ***************************************************************************
% Description: compute the cash flow under the situation of goodwill
% loss
%
% Parameters:
% T: planning horizon length
% d: (1*T) demands
% p: (1*T) prices
% s: (1*T) fixed ordering costs
% c: (1*T) unit vari ordering costs
% h: (1*T) unit holding costs
% B0: initial cash balance for the retailer
% beta: goodwill loss rate
% BL: credit-based loan
% TL: length of credit-based loan
% rL: interest rate of loan
%
%
% Decision variables:
% I: (1*T) end-of-period inventory of each period
% B: (1*T) end-of-period cash of each period
% BB: (T*T)end-of-period cash in each ordering cycle of the matrix
% WW: (T*T) end-of-period lost sale in each ordering cycle of the matrix
% x: (1*T) binary variables signaling whether ordering in each period
% y: (1*T) ordering quantity in each period
% w: (1*T) demand shortage (lost sales) in each period
%
% author: Zhen Chen
% time: 2019-02-19, 10:16
% ***************************************************************************

T = length(d);
BB = zeros(T, T); B = zeros(1, T); y = zeros(1, T);
B0 = B0 + BL;
loanPayBack = BL*(1 + TL)^rL; % quantity of loan and intersts needing to pay back
optIniB = zeros(1, T + 1); % optimal initial cash for each period
optIniB(1) = B0; 
optX = zeros(T, T);  % optimal ordering plan till each period t
satDemand = cell(T, T); % optimal satifying quantity in each optimal plan till period t 
XX = cell(T, T); % optimal ordering plan in an ordering round
WW = zeros(T, T); % lost sale quantity in the last period of an ordering cycle 
iniWIJ = zeros(T, T); % lost sale quantity in the last period of the first ordering cycle in an ordering round
iniBIJ = zeros(T, T); iniBIJ(1, :) = B0 * ones(1, T); % end-of-period cash in the last period of the first ordering cycle in an ordering round
whetherAdjustPlan = 0; whetherMoveOrderQuantity = 0; % recording whether ordering plan needs changing or inventory needs transfering


%% forward recursion
i = 1;
while i <= T
    for j = i : T
        % compute two ordering cycles in an ordering round
        if i > 1 && sum(optX(i - 1, : )) > 0
            for jj = i - 1 : -1 : 1
                if optX(i - 1, jj) == 1
                    m = jj; n = j; k = i - 1; % the ordering cycle: m ~ k, k + 1 ~ n
                    break;
                end
            end
            % record end-of-period cash lastB in the first ordering cycle, avoid cash decrease in the plan
            if i == j
                lastB = optIniB(k + 1); 
                iniSatDemand = satDemand{m, k}; lastSatDemand = iniSatDemand(end - k + m : end);
                lastX = zeros(1, k - m + 1); lastX(1) = 1; % lastX is the previous optimal ordering plan
            else
                lastB = BB(i, j - 1);
                lastSatDemand = satDemand{i, j - 1};
                lastX = XX{i, j - 1};
            end
            % judging whether needing to pay back the loan in this ordering
            % round
            if TL - m + 1 > 0
                thisL = L - m + 1;
                ratePay = loanPayBack;
            else
                thisL = 0; ratePay = 0;
            end
            xMn = zeros(1, n - m + 1); xMn(1) = 1; xMn(k + 2 - m) = 1; % two ordering cycles in an ordering round
            lastW = WW(i, j); % previous lost sale quantity in the last period of cycle i to j
            [~, xStart] = find(xMn == 1); xLast = diff(xStart); x1 = m + xStart(1) - 1; k1 = x1 + xLast - 1; % the first ordering cycle index: m, k1
            [iniB, iniW] = GetIniBW(xMn, iniBIJ(x1, k1), iniWIJ(x1, k1), BB(m : n, m : n), WW(m : n, m : n)); % get initial cash and initial lost sale in the ordering cycles of the ordering round
            [XX{i, j}, BB(i, j), satDemand{i, j}, finalIniB, finalIniW, WW(i, j)] = ComputeMkn(xMn, lastX, iniB, iniW, lastB, lastW, lastSatDemand, d(m : n), p(m : n), c(m : n), s(m : n),h(m : n), beta, 0, thisL, ratePay);
            if sum(XX{i,j}) > 1 % means new ordering cycle is added in the ordering round, so update initial cash and lost sale
                [iniBIJ(m : n, m : n), iniWIJ(m : n, m : n)] = UpdateIniBW(XX{i, j}, finalIniB, finalIniW, iniBIJ(m : n, m : n), iniWIJ(m : n, m : n));% 更新初始资金库存
            end
        end
        

    end
end

%% backward computing optimal ordering plan, satifying quantity and effective demand
x = optX(T,:); 
orderIndex = find(x == 1); 
satDemandFinal = zeros(1, T); 
xSum = sum(x);
k = T; % period index in the final period of the ordering round
if xSum > 0
    m = orderIndex(end); % period index in the last ordering cycle
    if sum(x(1 : m))==1
        satDemandFinal(m : k) = satDemand{m, k};
    end
    while sum(x(1 : m)) > 1 
        tempDemand = satDemand{m, k}; dLength = length(tempDemand);
        if dLength == k - m + 1
            satDemandFinal(m : k) = tempDemand;
            break;
        end
        % some times, tempDemand is longer than k - m + 1, because plan
        % adjustment
        if k == T
            satDemandFinal(k - dLength + 1 : k) = tempDemand;
        else
            satDemandFinal(k - dLength + 1 : m - 1) = tempDemand(1 : m + dLength - k - 1);
        end
        lastm = m;
        m= k - dLength + 1;
        k = lastm - 1;
    end
end
[Ed, ~, w] = ComputeEd(satDemandFinal, d, 0, beta);

%% compute y
i = 1;
xEndTo = zeros(1, T); % recording the index of the last period in an ordering cycle for each period
while i <= T
    if x(i) == 1 && i ~= T
        for j = i + 1 : T
            if x(j) ~= 1 && j < T
                continue
            else
                if x(j) == 1
                    endIndex = j - 1;
                end
                if j == T && x(j) ~= 1
                    endIndex = j;
                end
                y(i) = sum(satDemandFinal(i : endIndex));
                xEndTo(i : endIndex) = endIndex;
                i=endIndex;
                break;
            end
        end
    end
    if i==T
        if x(i) == 1
            y(T) = satDemandFinal(T);
            xEndTo(T) = T; 
        end
        break;
    end
    i = i + 1;
end

%% move order quantity
yIncrease = zeros(1, T);
for i = 1 : T - 1
    if x(i) == 1 
        if optIniB(i) - s(i) - c(i) * y(i) > 1e-4 && xEndTo(i) < T
            i2 = xEnd(i) + 1; % starting period of next ordering cycle
            if c(i) + sum(h(i : i2 - 1)) < c(i2)
                yIncrease(i) = min((optIniB(i2) - s(i2) - c(i2) * y(i2))/(c(i2) - sum(h(i : i2 - 1))),(optIniB(i) - s(i))/c(i) - y(i));
                y(i) = y(i) + yIncrease(i);
                y(i2) = y(i2) - yIncrease(i);
                optIniB(i2) = optIniB(i2) - yIncrease(i) * (c(i) + sum(h(i : i2 - 1)));
                whetherMoveOrderQuantity=1;
            end
        end   
    end
end

%% compute inventory flow and cash flow
I = zeros(1, T); I(1) = y(1) - (Ed(1) - w(1));
B(1) = B0 + p(1)*(Ed(1) - w(1)) - h(1) * I(1) - y(1) * c(1) - x(1)*s(1);
for t = 2 : T
    I(t) = I(t - 1) + y(t) - (Ed(t) - w(t));
    B(t) = B(t - 1) + p(t) * (Ed(t) - w(t)) - h(t) * I(t) - y(t) * c(t) - x(t) * s(t);
    if t == TL
        B(t) = B(t) - loanPayBack;
    end
end 
% B(T) = max(BB(1:T, T));
end