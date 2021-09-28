function [el,el_a]=exact_diag(n_body)
size_hilb_space= 2^n_body;
t_x = 2*ones(1,n_body);
lambda =0.2;
alpha =3;
sigma_x =[0, 1 ;...
          1, 0];
sigma_z = diag([1,-1]);
sigma_y =[0, -1i ;...
          1i, 0];
      
      
 name_long_range = 'thomas/a3_rcut30_p9_exp_param_n.mat';
 d_exp=load(name_long_range);
 
list_ham_xx = [];
list_ham_yy = [];
list_ham_zz = [];
i_h_nn=[];
[i_h, coup_power, coup_exp,coup_power_app,coup_nn] = long_range_coupling(n_body,alpha,lambda,d_exp);
for pos=1:n_body-1
i_h_nn{pos}=[pos, pos+1];
end
%neirest neighbours
%i_h=i_h_nn;
for k =1:length(i_h)
list_ham_xx{end+1} =reshape(kron(sigma_x, sigma_x) ,2,2,2,2);
list_ham_yy{end+1} =reshape(kron(sigma_y, sigma_y) ,2,2,2,2);
list_ham_zz{end+1} =reshape(kron(sigma_z, sigma_z) ,2,2,2,2);
end
OPTS.disp=1;
OPTS.isreal=0;
OPTS.issym=1;
[Ul,el,fl]= eigs(@L_v,size_hilb_space,4,'sr',OPTS);
diag(el)
%[Ul,el_a,fl]= eigs(@L_v_app,size_hilb_space,4,'sr',OPTS);
%diag(el_a)
function y=L_v(x)
        dim=size(x);
        
        x=reshape(x,t_x);
        r=zeros(t_x);
        %x=scon({x,G,A,A,cG},{[1,2,6,4],[1,-1,3],[3,-2,5,2],[5,-4,7,4],[6,-3,7]});
        c_v= -1 :(-1) : -n_body;
        for k=1:size(list_ham_xx,2)
            
            c_v_x=c_v;
            num_term=length(i_h{k});
            c_v_h = 1:num_term;
            for i_t =1 : num_term
                
                c_v_x(i_h{k}(i_t))=c_v_h(i_t);
                c_v_h(end+1)= -i_h{k}(i_t);
               
                
                
            end
             r=r+scon({lambda*coup_power(k)*x,list_ham_xx{k}},{c_v_x,c_v_h});
             r=r+scon({lambda*coup_power(k)*x,list_ham_yy{k}},{c_v_x,c_v_h});
             r=r+scon({lambda*coup_power(k)*x,list_ham_zz{k}},{c_v_x,c_v_h});
              r=r+scon({2*coup_exp(k)*x,list_ham_xx{k}},{c_v_x,c_v_h});
             r=r+scon({2*coup_exp(k)*x,list_ham_yy{k}},{c_v_x,c_v_h});
             r=r+scon({2*coup_exp(k)*x,list_ham_zz{k}},{c_v_x,c_v_h});

        end
        y=reshape(r,dim);
end
function y=L_v_app(x)
        dim=size(x);
        
        x=reshape(x,t_x);
        r=zeros(t_x);
        %x=scon({x,G,A,A,cG},{[1,2,6,4],[1,-1,3],[3,-2,5,2],[5,-4,7,4],[6,-3,7]});
        c_v= -1 :(-1) : -n_body;
        for k=1:size(list_ham_xx,2)
            
            c_v_x=c_v;
            num_term=length(i_h{k});
            c_v_h = 1:num_term;
            for i_t =1 : num_term
                
                c_v_x(i_h{k}(i_t))=c_v_h(i_t);
                c_v_h(end+1)= -i_h{k}(i_t);
               
                
                
            end
             r=r-scon({coup_power_app(k)*x,list_ham_xx{k}},{c_v_x,c_v_h});
             %r=r-scon({coup_power(k)*x,list_ham_yy{k}},{c_v_x,c_v_h});
        end
        y=reshape(r,dim);
end
    
end   
 
function [i_h, coup_power, coup_exp, coup_power_app,coup_nn] = ...
    long_range_coupling(n_body,alpha,lambda,d_exp)
i_h =[];
coup_power =[];
coup_power_app =[];
coup_exp =[];
coup_nn=[];
for position_1 = 1:n_body
    for position_2 = position_1+1: n_body
         distance = position_2-position_1;
         i_h{end+1}= [position_1, position_2];

         coup_power(end+1) = 1/distance^alpha;
         coup_power_app(end+1)=d_exp.Y(distance);
                  if distance >=2 
         coup_exp(end+1) = lambda^distance;
                  else
                      coup_exp(end+1) =0;
                      coup_nn(end+1)=lambda;
                  end
    end
end
end
