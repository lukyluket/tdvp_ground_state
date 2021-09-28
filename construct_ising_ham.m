function ham=construct_ising_ham(N,pos_def, theta)
sx=[0,1/2;1/2,0];
sz=diag([1/2,-1/2]);

bulk_ham=-cos(theta)*kron(sx,sx)-2*sin(theta)/2*(kron(sz,eye(2))+...
    kron(eye(2),sz)); 
% 
% def_ham =-sin(theta)/2*(kron(sz,eye(2))+...
%     kron(eye(2),sz)); 
% 
left_bound_ham = -2*sin(theta)/2*(kron(eye(2),sz));
rigth_bound_ham = -2*sin(theta)/2*(kron(sz,eye(2)));

for n=1:N-1
%     if n==pos_def
%         ham{n} = reshape(def_ham,2,2,2,2);    
%     else
    ham{n} = reshape(bulk_ham,2,2,2,2);    
    %end
end
ham{1} = ham{1} + reshape(left_bound_ham, 2,2,2,2);
ham{N-1} =ham{N-1} + reshape(rigth_bound_ham, 2,2,2,2);
end
