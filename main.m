function main

% *************************************************************************** 
% Description: the main function to compute the cash constrained lot sizing
% problem with credit-based loan and goodwill loss
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
% x: (1*T) binary variables signaling whether ordering in each period
% y: (1*T) ordering quantity in each period
% w: (1*T) demand shortage (lost sales) in each period
% 
% author: Zhen Chen
% time: 2019-02-14, 22:16
% ***************************************************************************


%% parameters
d = [203.18	134.21	77.91	200.79	 75.91	 51.12	225.93	321.30	182.27	171.49	128.64	101.63];
p = [16	 16  16	 17	 17	20	17	16	20	20	17	16];	
c = [13	13	13	 13	 13	 13	 13	 13	 13	 13	 13	 13];	
h = [1	1  1  1	 1	 1	  1	  1	  1	  1	  1	  1];	
s = [1000 1000 1000 1000 1000 1000 1000	1000 1000 1000 1000	1000];	
B0 = 6000;
beta = 0.1;
BL = 1000;
TL = 3;
rL = 0.05;
T = length(d);
pai = zeros(1, T); 

% %% or reading some parameter values from txt file in the current directory
% string = input('please input your file name: \n', 's'); % for example, d_12-84.txt
% %B0 = 2927; beta = 0.5;
% %BL = 2000; TL = 6; rL = 0.05;
% fidin = fopen(string);
% if fidin == -1
%     disp('没有这个数据文件\n');
% end
% index = 0;
% while ~feof(fidin)                        % 判断是否为文件末尾               
%     ch = fscanf(fidin, '%c', 1);  
%     if ch == ':'
%         index = index + 1;
%         switch index
%         case 1
%             T = fscanf(fidin, '%d', 1); p = zeros(1, T); c = zeros(1, T);    
%             h = zeros(1, T); s = zeros(1, T); d = zeros(1, T);     
%         case 2
%             for i = 1 : T
%                 p(i) = fscanf(fidin, '%f', 1);
%             end
%         case 3
%             for i = 1 : T
%                 c(i) = fscanf(fidin, '%f', 1);
%             end
%         case 4
%             for i = 1 : T
%                 h(i) = fscanf(fidin, '%f', 1);
%             end
%         case 5
%             for i = 1 : T
%                 s(i) = fscanf(fidin, '%f', 1);
%             end
%         case 6
%             for i = 1 : T
%                 d(i) = fscanf(fidin, '%f', 1);
%             end
%             break;
%         end  
%     end    
% end
% fclose(fidin);
% pai = zeros(1,T); 

%% run the forward algorithm
tic
whetherAdjustPlan = 0; whetherMoveOrderQuantity = 0;
if beta > 0
    [x, y, w, I, B, whetherAdjustPlan, whetherMoveOrderQuantity] = CashFlowGoodwill(d, p, s, c, h, beta, B0, BL, TL, rL); % there is no unit penalty cost for lost sale for this model
else
    [x, y, w, I, B] = CashFlowNoGoodwill(d, p, s, c, h, pai, B0, BL, TL, rL);
end
toc

%% output results to txt file
string = 'lot sizing-dp.txt';
OutputResult(x, y, w, I, B, beta, d, string, whetherAdjustPlan, whetherMoveOrderQuantity);

%% run the mip model by cplex
tic
fprintf ('\nRun MIP: \n');
[x, y, w, I, B] = MipCplex(d, p, s, c, h, pai, beta, B0, BL, TL, rL);
toc

%% output results to txt file
string = 'lot sizing-Cplex.txt';
OutputResult(x, y, w, I, B, beta, d, string, whetherAdjustPlan, whetherMoveOrderQuantity);

end