function [x, y, w, I, B, whetherAdjustPlan, whetherMoveOrderQuantity] = CashFlowGoodwill(d, p, s, c, h, beta, B0, BL, TL, rL, r0, Td)

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
% Td: payment delay length from customers
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
loanPayBack = BL*(1 + rL)^TL; % quantity of loan and intersts needing to pay back
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
        %% change [0 0 1] to [0 1 1]
        if i == 1 || sum(optX(i - 1, :)) == 0
            if i == 1   
                iniW = 0; % initial lost sale
            else
                [~, ~, tempW] = ComputeEd(zeros(1, i - 1), d(1 : i - 1), 0, beta);
                iniW = tempW(end); 
            end
            if i == j
                lastSatDemand = []; lastB = optIniB(i); lastX = [];
            else
                lastB = BB(i, j - 1); lastSatDemand = satDemand{i, j - 1}; lastX = XX{i, j - 1};
            end
            if sum(lastX) == 2 % compare with the last change && lastB > max(BB(1 : j - 1, j - 1)) - 1e-2 %上个阶段调整了,并且调整后是所在列最大, 要比较
                [~, xStart] = find(lastX == 1); xLength = length(lastX);
                m = j - xLength; n = j; k = m + xStart(2) - 1; % 最后一个 production cycle
                if TL - m + 1 > 0
                    thisL = TL - m + 1; ratePay = loanPayBack;
                else
                    thisL = 0; ratePay = loanPayBack;
                end
                iniB = zeros(1, 2); iniW = iniB;
                if m <= TL
                    iniB(1) = B0;
                else
                    iniB(1) = B0 - loanPayBack;
                end
                iniB(2) = iniBIJ(k, m + xLength - 1); iniW(1) = iniWIJ(m, k - 1); iniW(2) = iniWIJ(k, m + xLength - 1);
                xIJ = [lastX, 0]; preFinalW = [];
                [XX{i, j}, BB(i, j), satDemand{i, j}, finalIniB, finalIniW, WW(i, j)] = ComputeMkn(xIJ,...
                    lastX, iniB, iniW, lastB, preFinalW, lastSatDemand, d(m : n), p(m : n), c(m : n), s(m : n), h(m : n), beta, 0, thisL, loanPayBack);
                if i <= TL
                    iniB = B0;
                else
                    iniB = B0 - loanPayBack;
                end
                xIJ = zeros(1, j - i + 1); xIJ(1) = 1; iniW = tempW(end); preFinalW=[]; % 既定一个生产安排
                if TL - i + 1 > 0
                    thisL = TL - i + 1; ratePay = loanPayBack;
                else
                    thisL = 0; ratePay = 0;
                end
                [XX2, BB2, satDemand2, finalIniB2, finalIniW2, WW2] = ComputeMkn(xIJ, lastX, iniB, iniW, lastB, preFinalW,...
                                          lastSatDemand, d(i : j), p(i : j), c(i : j), s(i : j), h(i : j), beta, 0, thisL, ratePay);
                if BB2 > BB(i, j) + 1e-1    
                    XX{i, j} = XX2; BB(i, j) = BB2;
                    satDemand{i, j} = satDemand2;
                    finalIniB = finalIniB2; finalIniW = finalIniW2; WW(i, j) = WW2;
                end 
                tStart = j - length(XX{i, j}) + 1;
                [iniBIJ(tStart : n, tStart : n), iniWIJ(tStart : n, tStart : n)] = UpdateIniBW(XX{i, j},...
                           finalIniB, finalIniW, iniBIJ(tStart : n, tStart : n), iniWIJ(tStart : n, tStart : n)); % update initial lost sale and initial cash
            else
                if i <= TL
                    iniB = B0;
                else
                    iniB = B0 - loanPayBack;
                end
                xIJ = zeros(1, j - i + 1); xIJ(1) = 1; preFinalW = []; % set an ordering plan
                if TL - i + 1 > 0
                    thisL = TL - i + 1; ratePay = loanPayBack;
                else
                    thisL = 0; ratePay = 0;
                end
                [XX{i, j}, BB(i, j), satDemand{i, j}, iniBIJ(i, j), iniWIJ(i, j), WW(i, j)] = ComputeMkn(xIJ, lastX, iniB, iniW, lastB,...
                    preFinalW, lastSatDemand, d(i : j), p(i : j), c(i : j), s(i : j), h(i : j), beta, 0, thisL, ratePay);
            end 
        end
        
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
                thisLoanlength = TL - m + 1;
                ratePay = BL*(1 + rL)^TL;
            else
                thisLoanlength = 0; ratePay = 0;
            end
            xMn = zeros(1, n - m + 1); xMn(1) = 1; xMn(k + 2 - m) = 1; % two ordering cycles in an ordering round
            lastW = WW(i, j); % previous lost sale quantity in the last period of cycle i to j
            [~, xStart] = find(xMn == 1); orderLength = diff(xStart); k = m + xStart(1) - 1; k1 = k + orderLength - 1; % the first ordering cycle index: m, k1
            [iniB, iniW] = GetIniBW(xMn, iniBIJ(k, k1), iniWIJ(k, k1), BB(m : n, m : n), WW(m : n, m : n)); % get initial cash and initial lost sale in the ordering cycles of the ordering round
            [XX{i, j}, BB(i, j), satDemand{i, j}, finalIniB, finalIniW, WW(i, j)] = ComputeMkn(xMn, lastX, iniB, iniW, lastB, lastW, lastSatDemand, d(m : n), p(m : n), c(m : n),...
                s(m : n),h(m : n), beta, 0, thisLoanlength, ratePay);
            if sum(XX{i,j}) > 1 % means new ordering cycle is added in the ordering round, so update initial cash and lost sale
                [iniBIJ(m : n, m : n), iniWIJ(m : n, m : n)] = UpdateIniBW(XX{i, j}, finalIniB, finalIniW, iniBIJ(m : n, m : n), iniWIJ(m : n, m : n));% 更新初始资金库存
            end
        end
        
       %% heuristic adjustment, change [0, 0, 1] to [0 1 1]
        if i > 1 && sum(optX(i - 1, :)) == 0 && j - i < 4 % avoid too many
            iniX = XX{i, j};
            if sum(iniX) <= 1 
                iniSatDemand = satDemand{i, j};
                for m = max(1, i - 6) : i - 1  % means not necessary to change too many
                    k = i - 1; n = j;  % m~k, k + 1 ~ n,
                    if m == 1
                        thisIniW = 0;
                    else
                        [~, ~, tempW] = ComputeEd(zeros(1, m - 1), d(1 : m - 1), 0, beta);
                        thisIniW = tempW(end);
                    end
                    lastSatDemand = [zeros(k - m + 1, 1); iniSatDemand];
                    lastX = zeros(1, n - m + 1); lastX(i - m + 1 : n - m + 1) = iniX;
                    preFinalW = WW(i, j); lastB = BB(i,j);
                    xMn = zeros(1, n - m + 1); xMn(1) = 1; xMn(i - m + 1) = 1; 
                    if m <= TL
                        iniB = B0;
                    else
                        iniB = B0 - loanPayBack;
                    end
                    [iniB, iniW] = GetIniBW(xMn, iniB, thisIniW, BB(m : n, m : n), WW(m : n, m : n));
                    if TL - m + 1 > 0
                        thisL = TL - m + 1; ratePay = loanPayBack;
                    else
                        thisL = 0; ratePay = 0;
                    end
                    [tempXX, tempBB, tempSatDemand, finalIniB, finalIniW, WW(i,j)] = ComputeMkn(xMn, lastX, iniB, iniW, lastB,...
                                                     preFinalW, lastSatDemand, d(m : n), p(m : n), c(m : n), s(m : n), h(m : n), beta, 1, thisL, ratePay);
                    if tempBB > lastB + 1e-2
                        XX{i, j} = tempXX; BB(i, j) = tempBB; satDemand{i, j} = tempSatDemand;
                        [iniBIJ(m : n, m : n), iniWIJ(m : n, m : n)] = UpdateIniBW(XX{i, j}, finalIniB, finalIniW, iniBIJ(m : n, m : n), iniWIJ(m : n, m :n));
                    end
                end
            end
        end

    end % the loop for j
    
     %% compute optX
    if i > 1 && i < T
        [sortB , ia, ~] = unique(round(BB(1 : i, i))); % return no repetition of an array
        if length(sortB) > 1
            if (sortB(end) - sortB(end - 1)) / sortB(end) < 0.01 && WW(ia(end), i) > 10 && WW(ia(end -1 ), i) < 1 % if another one final cash very close, but lost sale is less, choose it
                startIndex = ia(end - 1); optIniB(i + 1) = BB(startIndex, i);
            else
                startIndex = ia(end); optIniB(i + 1) = BB(startIndex, i);
            end
        else
            [optIniB(i + 1), startIndex] = max(BB(1 : i, i)); 
        end
    else
        [optIniB(i + 1), startIndex] = max(BB(1 : i, i)); 
    end
    tempXX = XX{startIndex, i};
    if startIndex > 1
        optX(i, :) = optX(startIndex - 1, :);
    end
    optX(i, i - length(tempXX) + 1 : i) = tempXX; 
    
    %% heuristic adjustment, change plan [1 0 0 1 0] to [1 0 1 1 0]
    if i > 1       %启发式调整，%起始点 m~k,k+1~n  %调用子程序,情形1或许也能在这里调整
        n = i;
        if n == T % the last period does not need loss less
            needLossless = 0;
        else
            needLossless = 1;
        end
        if sum(optX(i, :)) > 1  % find two consecutive ordering cycle: m~k, k + 1~n 
            xNum = 0;
            for jj = i : -1 : 1
                if xNum == 0 && optX(i, jj) == 1
                    k = jj - 1; xNum = 1; continue;
                end
                if xNum == 1 && optX(i, jj) == 1
                    m = jj; break;
                end
            end
            % set the ordering plan
            xMk = zeros(1, k - m + 1); xMk(1) = 1; xKn = zeros(1, n - k); xKn(1) = 1; lastX = [xMk, xKn];
            lastB = BB(k+1, n); lastW = WW(k + 1, n); lastSatDemand = satDemand{k + 1, n};
            finalB = lastB; 
            % initial cash and lost sale
            [~, xStart] = find(lastX == 1); orderLength = diff(xStart); k = m + xStart(1) - 1; k1 = k + orderLength-1; % the first ordering cycle
            [iniB, iniW] = GetIniBW(lastX, iniBIJ(k, k1 ), iniWIJ(k, k1), BB(m : n, m : n), WW(m : n, m : n));
            for iChange = 2 : k - m + 1
                xMk = zeros(1, k - m + 1); xMk(1) = 1; xMk(1) = 1; xMk(iChange) = 1; xMn = [xMk, xKn]; 
                if TL - m + 1 > 0
                    thisLoanlength = TL - m + 1; loanPayBack = BL*(1 + TL)^rL;
                else
                    thisLoanlength = 0; loanPayBack = 0;
                end  
                [tempXX, tempBB, thisTempSatDemand, thisFinalIniB, thisFinalIniW, tempWW] = ComputeMkn(xMn, lastX, iniB, iniW, lastB,...
                    lastW, lastSatDemand, d(m : i), p(m : i), c(m : i), s(m : i), h(m : i), beta, needLossless, thisLoanlength, loanPayBack);
                if tempBB > finalB + 1e-2 % choose a bigger cash
                    BB(k + 1, n) = tempBB; WW(k + 1, n) = tempWW; XX{k + 1, n} = tempXX; optX(n, n - length(tempXX) + 1 : n) = tempXX;
                    finalB = tempBB; tempSatDemand = thisTempSatDemand; finalXX = tempXX;
                    finalIniB = thisFinalIniB; finalIniW = thisFinalIniW;
                end
            end % optX(k, m : k) = optX(i, m : k); i = k; % may use this 
            if finalB > lastB + 1e-4 % check whether [1 0 0 1 0 1 0 0 1 0] is less than  [1 0 0 0 1 0 1 0 0 0 0], sometimes not necessary
                satDemand{k + 1, n} = tempSatDemand;
                [iniBIJ(m : n, m : n), iniWIJ(m : n, m : n)] = UpdateIniBW(finalXX, finalIniB, finalIniW, iniBIJ(m : n, m : n), iniWIJ(m : n, m : n));
                for jj = n + 1 : T
                    xnj = zeros(1, jj - n); xMn = [XX{k + 1, n}, xnj];
                    [~, xStart] = find(xMn == 1); orderLength = diff(xStart); k = m + xStart(1) - 1; k1 = k + orderLength(1) - 1; % the first ordering cycle
                    [iniB, iniW] = GetIniBW(xMn, iniBIJ(k, k1), iniWIJ(k, k1), BB(m : jj, m : jj), WW(m : jj, m : jj));
                    lastSatDemand = satDemand{k + 1, jj}; lastX = XX{k + 1, jj}; lastW = WW(k + 1, jj); lastB = BB(k + 1, jj);
                    [tempXX, tempBB, thisTempSatDemand, thisFinalIniB, thisFinalIniW, tempWW] = ...
                        ComputeMkn(xMn, lastX, iniB, iniW, lastB, lastW, lastSatDemand, d(m : jj), p(m : jj), c(m : jj), s(m : jj),...
                                    h(m : jj), beta, 1, thisLoanlength, loanPayBack);
                    if tempBB > lastB + 1e-2 
                        BB(k + 1, jj) = tempBB; WW(k + 1, jj) = tempWW; XX{k + 1, jj} = tempXX; optX(k + 1, jj - length(tempXX) + 1 : jj) = tempXX;
                        satDemand{k + 1, jj} = thisTempSatDemand; finalXX = tempXX;
                        finalIniB = thisFinalIniB; finalIniW = thisFinalIniW;
                        [iniBIJ(m : jj, m : jj), iniWIJ(m : jj, m : jj)] = UpdateIniBW(finalXX, finalIniB, finalIniW, iniBIJ(m : jj, m : jj), iniWIJ(m : jj, m : jj));
                    end
                end
            end
        end
    end
    
    %% heuristic adjustment, change plan [1 0 0] to [0 1 0]
    if (sum(optX(i,:)) == 1) && i > 1
        [~, xStart] = find(optX(i, :) == 1);
        if xStart(1) == 1
            n = i;
            lastX = optX(i, 1 : i);
            lastB = BB(1,n); lastW = WW(1, n); lastSatDemand = satDemand{xStart(end), n};
            finalB = lastB;
            for iChange = xStart(1) + 1 : i
                m = iChange; xMn = zeros(1, n - m + 1); xMn(1) = 1;
                [~, ~, tempW] = ComputeEd(zeros(1, m - 1), d(1 : m - 1), 0, beta); iniW = tempW(end);
                if m <= TL
                    iniB = B0;
                else
                    iniB = B0 - loanPayBack;
                end
                if TL - m + 1 > 0
                    thisLoanlength = TL - m + 1;
                    ratePay = BL*(1 + TL)^rL;
                else
                    thisLoanlength = 0; ratePay = 0;
                end
                [thisTempXX, thisTempBB, thisTempSatDemand, thisFinalIniB, thisFinalIniW, thisTempWW] = ComputeMkn(xMn,...
                                         lastX, iniB, iniW, lastB, lastW, lastSatDemand, d(m : n), p(m : n), c(m : n), s(m : n), h(m : n), beta, 1, thisLoanlength, ratePay); 
                if thisTempBB > finalB + 1e-2 % choose a plan having more cash
                    tempXX = thisTempXX; finalB = thisTempBB; finalIniB = thisFinalIniB; finalIniW = thisFinalIniW;
                    tempSatDemand = thisTempSatDemand; tempWW = thisTempWW; iCh = iChange;
                end
            end
            if finalB > lastB + 1e-4
                BB(iCh,n) = finalB; XX{iCh,n} = tempXX; optX(n,:) = 0; optX(n, n - length(tempXX) + 1 : n) = tempXX;
                satDemand{iCh, n} = tempSatDemand; WW(iCh, n) = tempWW;
                [iniBIJ(m : n, m : n), iniWIJ(m : n, m : n)] = UpdateIniBW(xMn, finalIniB, finalIniW, iniBIJ(m : n, m : n),iniWIJ(m : n, m : n)); % update initial cash and initial lost sale
            end
        end
    end   
    i = i + 1; 
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
                i = endIndex;
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
            i2 = xEndTo(i) + 1; % starting period of next ordering cycle
            if c(i) + sum(h(i : i2 - 1)) < c(i2)
                yIncrease(i) = min((optIniB(i2) - s(i2) - c(i2) * y(i2))/(c(i2) - sum(h(i : i2 - 1))), (optIniB(i) - s(i))/c(i) - y(i));
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
r0Array = (1 / (1 + r0)) .^ (1 : T);
B(T) = B(T) / (1 + r0)^T;
% B(T) = max(BB(1:T, T));
end