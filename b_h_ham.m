function b_h_ham(d, mu,N)

a= create_a(d);
a_dagger = a';
num_op = a_dagger*a;

bulk_ham=-kron(a_dagger,a)+mu/2*(kron(num_op,eye(d))+...
    kron(eye(d),num_op));
left_bound_ham = mu/2*(kron(eye(d),num_op));
rigth_bound_ham = mu/2*(kron(num_op,eye(d)));
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


function a=create_a(d)
a=zeros(d);

for k=1:d-1
a(k,k+1)=sqrt(k);
end
end