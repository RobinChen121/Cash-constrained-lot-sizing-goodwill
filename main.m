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


d = [203.18	134.21	77.91	200.79	 75.91	 51.12	225.93	321.30	182.27	171.49	128.64	101.63];
p = [16	 16  16	 17	 17	20	17	16	20	20	17	16];	
c = [13	13	13	 13	 13	 13	 13	 13	 13	 13	 13	 13];	
h = [1	1  1  1	 1	 1	  1	  1	  1	  1	  1	  1];	
s = [1000 1000 1000 1000 1000 1000 1000	1000 1000 1000 1000	1000];	
B0 = 5000;
beta = 0.5;
BL = 1000;
TL = 6;
rL = 0.05;


I0 = 0; % initial inventory is zero
w0 = 0; % initial lost sales is zero


end