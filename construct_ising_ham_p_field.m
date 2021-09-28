function ham=construct_ising_ham_p_field(N, pos_def,c1,c2)
sx=[0,1;1,0];
sz=diag([1,-1]);

bulk_ham=kron(sx,sx)-c1/2*(kron(sz,eye(2))+...
    kron(eye(2),sz))-c2/2*(kron(sx,eye(2))+...
    kron(eye(2),sx)); 

def_ham =-c1/2*(kron(sz,eye(2))+...
    kron(eye(2),sz))-c2/2*(kron(sx,eye(2))+...
    kron(eye(2),sx)); 

left_bound_ham = -c1/2*(kron(eye(2),sz))-c2/2*(kron(eye(2),sx));
rigth_bound_ham = -c1/2*(kron(sz,eye(2)))-c2/2*(kron(sx,eye(2)));

for n=1:N-1
    if n==pos_def
        ham{n} = reshape(def_ham,2,2,2,2);    
    else
    ham{n} = reshape(bulk_ham,2,2,2,2);    
    end
end
ham{1} = ham{1} + reshape(left_bound_ham, 2,2,2,2);
ham{N-1} =ham{N-1} + reshape(rigth_bound_ham, 2,2,2,2);
end
