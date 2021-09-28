function[ l_c_mpo_left, l_c_mpo_rigth, l_O]= ...
    build_hamiltonian_xx(lambda,N)
      
[l_c_mpo_left{1}, l_c_mpo_rigth{1}, l_O{1}] = build_mpo(...
            [1,0,0,0;...
            0,0,lambda,0;...
            0,0,lambda,0;...
            0,0,0,1],...
            [0,-lambda,0,0;...
             0,0,0,0;...
             0,0,0,1;...
             0,0,0,0],...
            [0,0,0,0;...
             0,0,0,0;...
             0,0,0,0;...
             0,0,0,0],...
            [0,0,0,0;...
             0,0,0,0;...
             0,0,0,0;...
            0,0,0,0],...
           1,...
            N);
     
end
%        for k_exp = 2:length(d_exp.beta)
%              [l_c_mpo_left{k_exp}, l_c_mpo_rigth{k_exp}, l_O{k_exp}] = build_mpo(...
%              diag([1,d_exp.beta(k_exp),d_exp.beta(k_exp),d_exp.beta(k_exp),1]),...
%             [0,sin(theta)*d_exp.alpha(k_exp),0,0,0;...
%              0,0,0,0,1;...
%              0,0,0,0,0;...
%              0,0,0,0,0;...
%              0,0,0,0,0],...
%             [0,0,sin(theta)*d_exp.alpha(k_exp),0,0;...
%              0,0,0,0,0;...
%              0,0,0,0,1;...
%              0,0,0,0,0;...
%              0,0,0,0,0],...,...
%             [0,0,0,cos(theta)*d_exp.alpha(k_exp),0;...
%              0,0,0,0,0;...
%              0,0,0,0,0;...
%              0,0,0,0,1;...
%              0,0,0,0,0],...
%             1,...
%             N);
%         end
