function K = compute_kappa(A, rigth, H)
N = length(A);

K{N} = zeros(size(A{N-1},2),size(A{N-1},2)) ;
    
    for n = N-1: (-1) : 1
        K{n} = scon( {A{n},      conj(A{n}), K{n+1}},...
            {[-1, 1, 2],[-2, 3 ,2], [1,3]})+....
            scon({ A{n},       A{n+1},   conj(A{n}),  conj(A{n+1}),   H{n},    rigth{n+1} },...
            {[-1, 1, 3], [1, 2, 4], [-2, 7 , 8],   [7, 5, 6], [3,4,8,6], [2,5]});
    end

end