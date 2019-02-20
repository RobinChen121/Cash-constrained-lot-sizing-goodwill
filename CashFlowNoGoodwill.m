function [x, y, w, I, B] = CashFlowNoGoodwill(d, p, s, c, h, pai, B0, BL, TL, rL)

% ***************************************************************************
% Description: compute the cash flow under the situation of no goodwill
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
% pai: unit penalty cost for lost sale
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
% time: 2019-02-14, 22:16
% ***************************************************************************

T = length(d);
BB = zeros(T, T);
B = zeros(1, T);
x = zeros(1, T);
y = zeros(1, T);
w = zeros(1, T);
iniCash = B0 + BL;  % total initial cash at the beginnning of planning horizon
iniCashVector = zeros(1, T + 1); % maximum initial cash of each period (including period T + 1, which is the end-of-period cash of period T)
optCycleStartPeriod = zeros(1, T); % recording the optimal ordering cycle starting period in dynamic programming
options = optimset('TolFun', 1e-6, 'MaxIter', 50); % matlab function "linprog" setting (interior point algorithm by default)

%% forward dynamic programming 
for i = 1 : T
    for j = i : T
        cumUnitHoldCost = zeros(j- i + 1, j - i + 1); % cumulative unit holding cost for one ordering cycle
        tempH = zeros(j - i + 1, j - i + 1); % for computation of cumUnitHoldCost
        for m = 1 : j - i + 1
            for n = m : j - i + 1
                if n == m
                    tempH(m, n) = 0;
                else
                    tempH(m, n) = h( m + i - 1);
                end
            end
        end
        if j == i
            cumUnitHoldCost = 0;
        else
            cumUnitHoldCost(1,:) = tempH(1, :);
            for k = 2 : j - i + 1
                cumUnitHoldCost(k, :) = cumUnitHoldCost(k - 1, :) + tempH(k, :);
            end
        end
        f = c(i) * ones(1, j - i + 1) - p(i : j) - pai(i : j) + cumUnitHoldCost(j - i + 1, :); % coefficient vector of the objective function
        f=f';
        if i == 1
            iniCashVector(i) = iniCash;
        end
        lb = zeros( j - i + 1, 1); ub = d(i : j)'; % upper and lower bounds for decision variable (decision variable is the satisfying quantity in each period)
        A = zeros(j - i + 1, j - i + 1); b = zeros(j - i + 1, 1); % coefficient matrix of the constraints in the model
        b(j-i+1) = iniCashVector(i) - s(i);
        for k = 1 : j - i
            b(k) = iniCashVector(i) - s(i) - pai(i : k + i - 1) * d(i : k + i - 1)';
        end
        if TL - i + 1 > 0
            b(TL - i + 1 : j - i) = b(TL - i + 1 : j - i) - BL * (1 + rL)^TL;
        end
        A(j - i + 1, :)=c(i)*ones(1, j - i + 1);
        cumPrice = tril(repmat(p(i : j) + pai(i : j), n, 1));
        for k = 1 : j - i
            A(k, :) = c(i) * ones(1, j - i + 1) - cumPrice(k, :) + cumUnitHoldCost(k, :);
        end
        
        [~, fval, exitflag] = linprog(f, A, b, [], [], lb, ub, [], options); % use linprog to compute
        if exitflag == 1
            BB(i, j) = -fval - s(i) - pai(i:j) * d(i:j)' + iniCashVector(i);
            if BB(i,j) < iniCashVector(i) - pai(i:j)*d(i:j)' % if solution worse than all demandn lost, then no satifying quantity, demand all lost
                BB(i,j) = iniCashVector(i) - pai(i:j)*d(i:j)';
            end
        else
            BB(i, j) = iniCashVector(i) - pai(i : j) * d(i : j)'; % if no solution, then no satifying quantity, demand all lost
        end
        if i <= TL && j >= TL
            BB(i, j) = BB(i, j) - BL * (1 + rL)^TL;
        end
    end
    [iniCashVector(i + 1),optCycleStartPeriod(i)] = max(BB(1 : i, i)); % select an optimal ordering cycle
    % if cash increment is less than 0, optCycleStartPeriod(i) is
    % useless, make it negative
    if optCycleStartPeriod(i) <= TL && i >= TL
        if abs(iniCashVector(i + 1) + pai(optCycleStartPeriod(i) : i) * d(optCycleStartPeriod(i) : i)' + BL * (1 + rL)^TL - iniCashVector(optCycleStartPeriod(i))) < 1e-4
            optCycleStartPeriod(i) = -optCycleStartPeriod(i);
        end
    else
        if abs(iniCashVector(i + 1) + pai(optCycleStartPeriod(i) : i) * d(optCycleStartPeriod(i) : i)' - iniCashVector(optCycleStartPeriod(i))) < 1e-4
            optCycleStartPeriod(i) = -optCycleStartPeriod(i);
        end
    end
end

%% compute optimal satisfying quantity in the optimal plan
% backward searching to get the optimal plan x
j = T;
markUnfeasibleStartPeriod = zeros(1, T); % record the unfeasible cycle starting period, binary
while j >= 1
    if optCycleStartPeriod(j) > 0
        x(optCycleStartPeriod(j)) = 1;
        j = optCycleStartPeriod(j);
    else
        markUnfeasibleStartPeriod(-optCycleStartPeriod(j)) = 1;
        j = -optCycleStartPeriod(j);
    end
    j = j - 1;
end

% forward compute satisfying quantity
satDemands = zeros(1, T);
i = 1;
while i <= T
    if x(i) == 1 && i ~= T
        % find the ordering cycle ending period
        for j = i + 1 : T
            if x(j) ~= 1 && j < T && markUnfeasibleStartPeriod(j) ~= 1
                continue;
            else
                if x(j) == 1
                    cycleEndPeriod = j - 1;
                end
                if markUnfeasibleStartPeriod(j) == 1
                    cycleEndPeriod = j - 1;
                end
                if j == T && x(j) ~= 1 && markUnfeasibleStartPeriod(j) ~= 1
                    cycleEndPeriod = T;
                end
                
                cumUnitHoldCost = zeros(cycleEndPeriod - i + 1, cycleEndPeriod - i + 1);
                tempH = zeros(cycleEndPeriod - i + 1, cycleEndPeriod - i + 1);
                for m = 1 : cycleEndPeriod - i + 1
                    for n = m : cycleEndPeriod - i + 1
                        if n == m
                            tempH(m, n) = 0;
                        else
                            tempH(m, n) = h(m + i - 1);
                        end
                    end
                end
                if cycleEndPeriod == i
                    cumUnitHoldCost = 0;
                else
                    cumUnitHoldCost(1, :) = tempH(1, :);
                    for k = 2 : cycleEndPeriod - i + 1
                        cumUnitHoldCost(k, :) = cumUnitHoldCost(k - 1, :)+tempH(k, :);
                    end
                end
                f = c(i) * ones(1,cycleEndPeriod - i + 1) - p(i : cycleEndPeriod) - pai(i : cycleEndPeriod) + cumUnitHoldCost(cycleEndPeriod - i + 1, :);
                f = f';
                lb = zeros(cycleEndPeriod - i + 1, 1); ub = d(i : cycleEndPeriod)';
                A = zeros(cycleEndPeriod - i + 2, cycleEndPeriod - i + 1); b = zeros(cycleEndPeriod - i + 2, 1);
                b(cycleEndPeriod - i + 2) = iniCashVector(i) - s(i);
                for k = 1 : cycleEndPeriod - i + 1
                    b(k) = iniCashVector(i) - s(i) - pai(i : k + i - 1) * d(i : k + i - 1)';
                end
                if TL - i + 1 > 0
                    b(TL - i + 1 : cycleEndPeriod - i) = b(TL - i + 1 : cycleEndPeriod - i) - BL*(1+rL)^TL;
                end
                A(cycleEndPeriod - i + 2, :) = c(i) * ones(1, cycleEndPeriod - i + 1);
                cumPrice = tril(repmat(p(i : cycleEndPeriod) + pai(i : cycleEndPeriod), n, 1));
                for k = 1 : cycleEndPeriod - i + 1
                    A(k, :) = c(i) * ones(1, cycleEndPeriod - i + 1) - cumPrice(k, :) + cumUnitHoldCost(k, :);
                end
                [satDemands(i : cycleEndPeriod), ~ , ~] = linprog(f, A, b, [], [], lb, ub, [], options);
                i=j;
                break;
            end
        end
    end
    if i == T
        if x(i) == 1
            satDemands(T) = min((iniCashVector(i) - s(i))/c(T), d(T)); % satisfying quantity in this situation can be computated directly
        end
        break;
    end
    if x(i) == 0
        i = i + 1;
    end
end

% get w
for i = 1 : T
    w(i) = d(i) - satDemands(i);
end


% get y
i = 1;
while i <= T
    if x(i) == 1 && i ~= T
        for j = i + 1 : T
            if x(j) ~= 1 && j < T
                continue
            else
                if x(j) == 1
                    cycleEndPeriod = j - 1;
                end
                if j == T && x(j) ~= 1
                    cycleEndPeriod = j;
                end
                y(i) = sum(satDemands(i : cycleEndPeriod));
                i = j;
                break;
            end
        end
    end
    if i == T
        if x(i) == 1
            y(T) = satDemands(T);
        end
        break;
    end
    if x(i) == 0
        i = i + 1;
    end
end

%% move ordering quantity
cycleStartPeriod = zeros(1, T);
cycleEndPeriod = zeros(1, T);
for i = 1 : T
    if x(i) == 1
        cycleStartPeriod(i) = i;
    end
end
for i = 1 : T
    if x(i) == 1
        for j = i + 1 : T
            if x(j) == 1
                cycleEndPeriod(i) = j - 1;
                break;
            end
            if j == T
                cycleEndPeriod(i) = j;
            end
        end
        if i == T
            cycleEndPeriod(i) = T;
        end
    end
end

yIncrease = zeros(1, T);
for i = 1 : T - 1
    if x(i) == 1 
        if iniCashVector(i) - s(i)- c(i) * y(i) > 1e-4 && cycleEndPeriod(i) < T && ...
                iniCashVector(cycleEndPeriod(i) + 1) - s(cycleEndPeriod(i) + 1) - c(cycleEndPeriod(i) + 1) * y(cycleEndPeriod(i) + 1) > 1e-4
            if c(i) + sum(h(i : cycleEndPeriod(i))) < c(cycleEndPeriod(i) + 1)
                yIncrease(i) = (iniCashVector(i) - s(i))/c(i) - y(i);
                y(i) = y(i) + yIncrease(i);
                y(cycleEndPeriod(i) + 1) = y(cycleEndPeriod(i) + 1) - yIncrease(i);
                iniCashVector(cycleEndPeriod(i) + 1) = iniCashVector(cycleEndPeriod(i) + 1) - yIncrease(i) * (c(i) + sum(h(i : cycleEndPeriod(i))));
                for j = cycleEndPeriod(cycleEndPeriod(i) + 1) + 1 : T + 1
                    iniCashVector(j) = iniCashVector(j) + yIncrease(i) * (c(cycleEndPeriod(i) + 1) - c(i) - sum(h(i : cycleEndPeriod(i))));
                end
            end
        end
    end
end

%% compute inventory flow and cash flow
I = zeros(1, T);
I0 = 0;  % initial inventory is zero
I(1) = I0 + y(1) - (d(1) - w(1));
if I(1) < 1e-3  % remove the decimal part
    I(1) = 0;
end
for t = 2 : T
    I(t) = I(t - 1) + y(t) - (d(t) - w(t));
    if I(t) < 1e-3
        I(t) = 0;
    end
end
B(1) = iniCash + p(1) * (d(1) - w(1)) - h(1) * I(1) - y(1) * c(1) - x(1) * s(1) - pai(1)*w(1);
for t = 2 : T
    if t == TL  % need to pay back the loan and interest at this period
        B(t) = B(t - 1) + p(t) * (d(t) - w(t))- h(t) * I(t) - y(t) * c(t) - x(t) * s(t) - pai(t) * w(t) - BL * (1 + rL)^TL;
    else
        B(t) = B(t -1 ) + p(t) * (d(t) - w(t))- h(t) * I(t) - y(t) * c(t) - x(t) * s(t) - pai(t)* w(t);
    end
end

end

