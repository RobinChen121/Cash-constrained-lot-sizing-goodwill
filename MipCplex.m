function [x, y, w, I, B] = MipCplex(d, p, s, c, h, ~, beta, B0, BL, TL, rL, r0, Td)

% ***************************************************************************
% Description: compute the cash flow via Cplex
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
% pai: unit penalty cost for lost sale, in the mip model, there is no pai
% r0: deposite interest rate
% Td: payment delay length from customers

% author: Zhen Chen
% time: 2019-02-18, 15:16
% ***************************************************************************

T = length(d);
%% compute the matrix values for the constraints
cumH = zeros(1, T);
cumD = zeros(1, T);
for i = 1 : T
    cumH(i) = sum(h(i : T));
end
for i = 1 : T
    cumD(i) = sum(d(1 : i));
end
cycleCumH = zeros(T, T);
for i = 1 : T
    for j = 1 : T
        if j >= i
            cycleCumH(i, j) = 0;
        else
            cycleCumH(i, j) = sum(h(j : i-1));
        end
    end
end
cyclePrice = repmat(p, T, 1);
cyclePrice = tril(cyclePrice, -Td - 1);
cycleBeta = zeros(T, T);
for i = 1 : T
    for j = 1 : T
        if i > 1 && i - j == 1
            cycleBeta(i, j) = beta;
        end
    end
end
tempOne=diag(ones(1, T));
tempOne(1,1) = 0;


%% f, Aineq, bineq
% decision variables: x, y , w, Ed, delta
r0Array = (1 / (1 + r0)) .^ (1 : Td);
r0Array2 = [ones(1, T - Td), r0Array];
f = [s , (cumH + c) , (cumH + p .* r0Array2), (-cumH - p .* r0Array2), zeros(1, T)]' / (1 + r0)^T; % objective function
M1 = ones(1, T) * sum(d);
M2 = sum(d);
sumB0 = BL + B0;
Aineq = [diag(-M1), diag(ones(1, T)), zeros(T, T), zeros(T, T),zeros(T, T); % yt <= Mxt
    tril(repmat(s, T, 1)), tril(repmat(c, T, 1)) + tril(cycleCumH), cyclePrice + tril(cycleCumH), -cyclePrice - tril(cycleCumH), zeros(T,T); % Ctyt + Stxt <= Bt-1 
    %zeros(T, T), diag(ones(1, T)), zeros(T, 3 * T); % traditional capacity model
    zeros(T, T), -tril(ones(T, T)), -tril(ones(T, T)), tril(ones(T, T)), zeros(T, T); % It>=0
    zeros(T, 2 * T), diag(ones(1, T)), -diag(ones(1, T)), zeros(T, T); % wt <= Edt
    zeros(T, 2 * T), cycleBeta, diag(ones(1, T)), M2 * tempOne;    %Edt �� dt - �� * w(t-1) + M2*(1-��t)
    zeros(T, 2 * T), -cycleBeta, -diag(ones(1, T)), +M2 * tempOne;    %Edt >= dt - �� * w(t-1) + M2 * (1-��t) - sigma
    zeros(1, 4 * T), [1, zeros(1, T - 1)] %��_1=1
    zeros(1, 4 * T),[-1, zeros(1, T - 1)] %��_1=1
    zeros(T, 3 * T), diag(ones(1, T)), -diag(d);   % Ed_t �� dt *��_t
    zeros(T - 1, 2 * T), -cycleBeta(2 : T, :),zeros(T - 1, T), -M2 * tempOne(2 : T, :);    %d_t �� ��w(t-1) + M2 *��_t
    zeros(T - 1, 2 * T), cycleBeta(2 : T, :), zeros(T - 1, T), M2 * tempOne(2 : T, :);     %d_t >=��w(t-1) + M2 * (1-��_t)
    ];
bineq = [zeros(T, 1); [sumB0 * ones(TL, 1); (sumB0 - BL * (1 + rL)^TL) * ones(T - TL, 1)]; zeros(T, 1); zeros(T, 1); d' + [0; M2 * ones(T - 1, 1)];
   -d' + [0; M2 * ones(T - 1, 1)]; 1; -1; zeros(T, 1); -d(2 : T)'; d(2 : T)' + M2 * ones(T - 1, 1)];

%    bineq = [zeros(T, 1);  zeros(T, 1); zeros(T, 1); d' + [0; M2 * ones(T  -1, 1)];   
%    -d' + [0; M2 * ones(T - 1, 1)]; 1; -1; zeros(T, 1); -d(2 : T)'; d(2 :  T)' + M2 * ones(T - 1, 1)];  % Aksen's model

%    bineq = [zeros(T, 1);  90 * ones(T, 1); zeros(T, 1); zeros(T, 1); d' + [0; M2 * ones(T  -1, 1)];   % traditional capacity model
%    -d' + [0; M2 * ones(T - 1, 1)]; 1; -1; zeros(T, 1); -d(2 : T)'; d(2 :  T)' + M2 * ones(T - 1, 1)];

Aeq = [];
beq = [];
lb = zeros(5 * T,1);
ub = [ones(T, 1); M2 * ones(2 * T, 1); d'; ones(T, 1)];

%ub = [ones(T, 1); M2 * ones(1 * T, 1); zeros(T, 1); d'; ones(T, 1)]; % traditional capacity model

% type of decision variables
a = zeros(1, 5 * T);
ctype = char(a);
for i = 1 : T
    ctype(i) = 'I';
end
for i = T + 1 : 2 * T
    ctype(i) = 'C';
end
for i = 2 * T + 1 : 3 * T
    ctype(i) = 'C';
end
for i = 3 * T + 1 : 4 * T
    ctype(i) = 'C';
end
for i = 4 * T + 1 : 5 * T
    ctype(i) = 'I';
end

%% solve
options = cplexoptimset;
options.Display = 'off';
options.MaxIter = 750000;
options.MaxTime = 18000;
[x, fval, ~, ~] = cplexmilp(f, Aineq, bineq, Aeq, beq,...
    [ ], [ ], [ ], lb, ub, ctype, [ ], options);
%fprintf ('Solution value = %f \n', -fval + (B0 + BL - BL * ((1 + rL)^TL))/ (1 + r0)^T);

%% output x, y, w, I, B
for i = 1 : T
    if abs(x(i)) < 1e-4
        x(i) = 0;
    end
end
Ed = zeros(1, T); % effective demands
Ed(1) = d(1);
for i = 2 : T
    Ed(i) = max(0, d(i) - beta * x(2 * T + i - 1));
end
y = zeros(1, T);
for i = T + 1 : 2 * T
    y(i - T) = x(i);
end
w = zeros(1, T);
for i = 2 * T + 1 : 3 * T
    w(i - 2 * T) = x(i);
end
% inventory flow
I = zeros(1, T);
I(1) = y(1) - (Ed(1) - w(1));
for t = 2 : T
    I(t) = I(t - 1) + y(t) - (Ed(t) - w(t));
end

B = zeros(1, T);

if Td == 0
    B(1) = sumB0 + p(1) * (Ed(1) - w(1)) - h(1) * I(1) - y(1) * c(1) - x(1) * s(1);  
else
    B(1) = sumB0 - h(1) * I(1) - y(1) * c(1) - x(1) * s(1);
end
for t = 2 : T
    if t <= Td
        B(t) = B(t - 1) - h(t) * I(t) - y(t) * c(t) - x(t) * s(t);
    else
        B(t) = B(t - 1) + p(t - Td) * (Ed(t - Td) - w(t - Td)) - h(t) * I(t) - y(t) * c(t) - x(t) * s(t);
    end
    if t == TL
        B(t) = B(t) - BL * (1 + rL)^TL;
    end
    if t == T
        delayPayPre = p(T - Td + 1 : T) .* (Ed(T - Td + 1 : T) - w(T - Td + 1 : T));
        delayPayFinal = delayPayPre .* (1 / (1 + r0)).^(1 : Td);
        delayPayFinal = sum(delayPayFinal);
        B(t) = B(t) + delayPayFinal;
    end
end
B(T) = B(T) / (1 + r0)^T;
x = x(1 : T);
end