function left_mpo = compute_left_mpo( A, O, cl_O_left,left)
N = length(A);
left_mpo{1} = cl_O_left;
left_mpo{2} = scon({left{1},A{1},    conj(A{1}), O{1},         cl_O_left},...
    {[2,5],[2,-1,3],[5,-3,4] , [1, -2, 3, 4], [-4,1]});
for n = 3: N
    left_mpo{n} = scon({A{n-1},    conj(A{n-1}), O{n-1},        left_mpo{n-1}},...
        {[4,-1,5],[1,-3,3],[2,-2,5,3],[4,2,1]});
end
end
