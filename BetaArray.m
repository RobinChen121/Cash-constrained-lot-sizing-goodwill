function a = BetaArray(n, beta)

% ***************************************************************************
% Description: get an array of beta for use in the constraint matrix
%
% Parameters:
% n: length of the array
% beta: goodwill loss rate
%
% Decision variables:
% a: (1*n) an array of beta
%
% author: Zhen Chen
% time: 2019-02-20, 09:16
% ***************************************************************************

a = zeros(1,n);
a(1) = beta^(n-1);
for i = 2 : n
    a(i) = a(i - 1)/beta;
end
for i = n : -1 : 1
    a(i) = a(i) * (-1)^(n - i);
end


end