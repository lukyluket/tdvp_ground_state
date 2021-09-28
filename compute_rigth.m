function rigth = compute_rigth(A,rigth)
L =  length(A);
%rigth{L} = eye(size(A{L},2))/size(A{L},2);
for k =L-1:-1:1
    
    rigth{k} = scon( {  A{k+1},    conj(A{k+1}),   rigth{k+1}},...
        {[-1, 1, 2],   [-2 , 3 , 2],     [1, 3]});
end

end
