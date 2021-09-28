function rigth_mpo = compute_rigth_mpo_l( A, O, cl_O_rigth,rigth_mpo,n)
N = length(A);
if n ==N
rigth_mpo{N} = cl_O_rigth;
elseif n==N-1
rigth_mpo{N-1} = scon({A{N},    conj(A{N}), O{N},         cl_O_rigth},...
    {[-1,2,3],[-3,2,4] , [-2, 1, 3, 4], [1,-4]});
%for n = N-2:(-1): 1
else
rigth_mpo{n} = scon({A{n+1},    conj(A{n+1}), O{n+1},        rigth_mpo{n+1}},...
        {[-1,4,5],[-3,1,3],[-2,2,5,3],[4,2,1]});
end
end