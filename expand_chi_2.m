function NA=expand_chi_2(A,O,chi_max,left,rigth,left_mpo, rigth_mpo,c_mpo_rigth, epsilon,delta_t)
N= length(A);

for n=N-1:(-1):2
    Vtr = fix_gauge_for_B_l(A, rigth,n);
    Vtl = fix_gauge_for_B_left_l(A, left,n-1);
    %Vtrr = fix_gauge_for_B_l (A, rigth,n);
   % A{n} = update_A_mpo_l(A, O,left, rigth, left_mpo, rigth_mpo, Vtrr,epsilon,delta_t,n);
    sileft = diag(sqrt(diag(pinv(left{n-1}, epsilon))));
    sirigth =diag(sqrt(diag(pinv(rigth{n}, epsilon))));
    tn={...
        sileft,...
        left_mpo{n-1},...
        A{n-1},...
        O{n-1},...
        O{n},...
        A{n},...
        rigth_mpo{n},...
        sirigth,...
        Vtr{n},...
        Vtl{n-1}};
    
    c_net={...
    [1,14],...%sileft
    [2,3,1],...%left_mpo{n-1}
    [2,5,4],...%A{n-1}
    [3,7,4,13],...%O{n-1}
    [7,9,6,12],...%O{n}
    [5,8,6],...%A{n}
    [8,9,10],...%rigth_mpo{n}
    [10,11],...%sirigth
    [-2,11,12],...% Vtr{n}
    [14,-1,13]...%Vtl{n-1}
    };


    xy = scon (tn,c_net);
    
    [u,s,v]= svd(xy);
    
    chi_l = size(A{n},1);
    chi_r = size(A{n},2);
    d= size(A{n},3);
    real_length = min(n-1,N-n+1);
    real_size =d^real_length;
    if chi_l+chi_max > real_size
    else    
    chi_t=min(real_size, chi_max);
    
    x= u*s(:,1:chi_t);
    vd=v';
    y= vd(1:chi_t,:);
    
    
    
    n_chi_l= size(y,1);
    n_chi_r = size(y,2);
    A{n} = update_A_mpo_l(A, O,left, rigth, left_mpo, rigth_mpo, Vtr,epsilon,delta_t,n);
    A{n}(chi_l+1:chi_l+n_chi_l, 1:chi_r,:) = scon({ y, conj(Vtr{n}), sirigth},...
            {[-1,1], [1,2,-3],[-2,2]})*sqrt(delta_t);
    
    chi_l = size(A{n-1},1);
    chi_r = size(A{n-1},2);
    d= size(A{n-1},3);
    
    n_chi_l= size(y,1);
    n_chi_r = size(y,2);
    A{n-1}(1:chi_l,chi_r+1:chi_r+n_chi_l,:) = scon({ sileft, conj(Vtl{n-1}), x},...
            {[-1,1], [1,2,-3],[2,-2]})*sqrt(delta_t);
        
      [A,left,rigth] = fix_rigth_to_identity_l (A, left, rigth, epsilon,n-1);
      rigth_mpo= compute_rigth_mpo_l(A, O, c_mpo_rigth,rigth_mpo,n-1);
    end
    end
NA=A;
end