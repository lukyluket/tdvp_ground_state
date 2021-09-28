function left = compute_left(A,left)
L =  length(A);
%left{1} = eye(size(A{1},1));
for k = 2: L+1
left{k} = scon( {  A{k-1},    conj(A{k-1}),   left{k-1}},...
        {[1, -1, 2],   [3 , -2 , 2],     [1, 3]});
end
end