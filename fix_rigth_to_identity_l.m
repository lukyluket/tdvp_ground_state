function [A, left, rigth, rigth_change] = fix_rigth_to_identity_l(A, left, rigth,epsilon,n)
%N = length(A);
%rigth{N} = eye( size(A{N},2));
if n >0
rigth{n} =   scon( {  A{n+1},    conj(A{n+1}),   rigth{n+1}},...
        {[-1, 1, 2],   [-2 , 3 , 2],     [1, 3]});
    
    [Ur, er] = eig(0.5*(rigth{n}+rigth{n}'));
    left_change  = sqrt(pinv(er,epsilon))*Ur';
    rigth_change = Ur*sqrt(er);
    A{n+1} = scon({A{n+1}, left_change},{[1,-2,-3],[-1,1]});
    A{n} = scon({A{n}, rigth_change},{[-1,1,-3],[1,-2]});
rigth{n} = scon( {  A{n+1},    conj(A{n+1}),   rigth{n+1}},...
        {[-1, 1, 2],   [-2 , 3 , 2],     [1, 3]});
    new_rigth = scon( {  A{n},    conj(A{n}),   rigth{n}},...
        {[-1, 1, 2],   [-2 , 3 , 2],     [1, 3]});
if n >1
    err_gf = abs(new_rigth - rigth{n-1});
    max_err = max(err_gf(1:end)); 
end
    
else
last_rigth =   scon( {  A{1},    conj(A{1}),   rigth{1}},...
        {[-1, 1, 2],   [-2 , 3 , 2],     [1, 3]});
    [Ur, er] = eig(0.5*(last_rigth+last_rigth'));
    left_change  = sqrt(pinv(er,epsilon))*Ur';
    rigth_change = Ur*sqrt(er);
    
    A{1} = scon({A{1}, left_change},{[1,-2,-3],[-1,1]});
    left{1} = scon({rigth_change, left{1}, conj(rigth_change)},{[1,-1],[1,2],[2,-2]});
end
end