function OutputResult(x, y, w, I, B, beta, d, string, whetherAdjustPlan, whetherMoveOrderQuantity)

% ***************************************************************************
% output result to a txt file

% author: Zhen Chen
% time: 2019-02-14, 22:16
% ***************************************************************************

fid = fopen(string, 'w');
fprintf(fid, 'whether ordering in this period: \n'); fprintf(fid,'x = ');
T = length(y);
for i = 1 : T
    fprintf(fid, '%d\t', x(i));
    if mod(i, 20)==0
        fprintf(fid, '\n'); % new line when planning length is very long
    end
end
fprintf(fid, '\n'); fprintf(fid, 'ordering quantity: \n'); fprintf(fid, 'y = ');
for i=1:T
    fprintf(fid, '%6.2f\t', y(i));
    if mod(i, 20) == 0
        fprintf(fid, '\n');
    end
end
fprintf(fid, '\n'); fprintf(fid, 'lost sales: \n'); fprintf(fid, 'w = ');
for i = 1 : T
    fprintf(fid, '%6.2f\t', w(i));
    if mod(i, 20) == 0
        fprintf(fid, '\n');
    end
end
fprintf(fid,'\n');fprintf(fid,'end-of-period inventory:\n'); fprintf(fid, 'I = ');
for i = 1 : T
    fprintf(fid, '%6.2f\t', I(i));
    if mod(i, 20) == 0
        fprintf(fid, '\n');
    end
end
fprintf(fid,'\n');fprintf(fid,'demand quantity in each period:\n');fprintf(fid, 'd = ');
for i = 1 : T
    fprintf(fid,'%6.2f\t', d(i));
    if mod(i,20) == 0
        fprintf(fid, '\n');
    end
end
Ed = zeros(1, T); Ed(1) = d(1); 
for i = 2 : T
    Ed(i) = max(0, d(i) - beta*w(i - 1)); % effective demand 
end
fprintf(fid, '\n'); fprintf(fid,'effective demand in each period:\n'); fprintf(fid, 'Ed = ');
for i = 1 : T
    fprintf(fid, '%6.2f\t', Ed(i));
    if mod(i, 20) == 0
        fprintf(fid, '\n');
    end
end
fprintf(fid, '\n'); fprintf(fid, 'satisfying demand in each period:\n'); fprintf(fid, 'SatD = ');
for i = 1 : T
    fprintf(fid, '%6.2f\t', Ed(i) - w(i));
    if mod(i,20) == 0
        fprintf(fid, '\n');
    end
end
fprintf(fid, '\n'); fprintf(fid, 'end-of-period cash in each period:\n'); fprintf(fid, 'B = ');
for i = 1 : T
    fprintf(fid, '%6.2f\t', B(i));
    if mod(i, 20) == 0
        fprintf(fid, '\n');
    end
end
fprintf ('Solution value of the forward algorithm is = %f \n', B(T));
if beta>0
    if whetherAdjustPlan == 1
        fprintf ('adjust the plan\n');
    end
    if whetherMoveOrderQuantity == 1
        fprintf ('move the ordering quantity\n');
    end
end
fclose(fid); 
end