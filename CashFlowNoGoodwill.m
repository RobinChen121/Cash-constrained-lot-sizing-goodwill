function [x,y,w,I,B] = CashFlowNoGoodwill(d, p, s, c, h, pai, B0, BL, TL, rL)

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
% BB: (T*T) maximum cash increment from period m to period n
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

for i = 1 : T
    for j = 1 : T
        
        
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
       B(t) = B(t - 1) + p(t) * (d(t) - w(t))- h(t) * I(t) - y(t) * c(t) - x(t) * s(t) - pai(t) * w(t) - BL * (1 + TL)^rL; 
    else
        B(t) = B(t -1 ) + p(t) * (d(t) - w(t))- h(t) * I(t) - y(t) * c(t) - x(t) * s(t) - pai(t)* w(t);
    end
end  


end

