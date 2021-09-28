function[ l_c_mpo_left, l_c_mpo_rigth, l_O]= ...
    build_hamiltonian_xx_yy_zz_exp_2(lambda,N,coup_const)
 last=8;     
operator_id =zeros(last);
operator_id(1,1) = 1;
operator_id(2,3) = lambda;
operator_id(3,3) = lambda;
operator_id(4,5)=  lambda;
operator_id(5,5) = lambda;
operator_id(6,7) = lambda;
operator_id(7,7) = lambda;
operator_id(last,last) = 1;

operator_sx= zeros(last);
operator_sx(1,2) = lambda;
operator_sx(3,last) = 1;

operator_sy =zeros(last);
operator_sy(1,4) = lambda;
operator_sy(5,last) = 1;

operator_sz =zeros(last);
operator_sz(1,6) = lambda;
operator_sz(7,end) =1;


[l_c_mpo_left{1}, l_c_mpo_rigth{1}, l_O{1}] = build_mpo(...
           operator_id,...
           operator_sx,...
           operator_sy,...
           operator_sz,...
           coup_const,...
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
