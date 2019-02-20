function [finalXMn, finalB, satDemand, finalIniB, finalIniW, finalW] = ComputeMkn(xMn, lastX, thisIniB, ... 
                                    thisIniW, lastB, lastW, lastSatDemand, d, p, c, s, h, beta, needLossless, thisLoanLength, ratePay)
                                
% ***************************************************************************
% Description: get objective matrix, constraint matrix for an ordering
% round to compute the maximum cash increment
%
% Parameters:
% xMn: ordering plan in an ordering cycle from period m to period n
% lastX: an array, previous optimal ordering plan
% thisIniB: initial cash for this ordering round
% thisIniW: initial lost sale for this ordering round
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
                                
                                
                                

end