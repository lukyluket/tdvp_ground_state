function NA = update_A_mpo_l(A, O,left, rigth, left_mpo, rigth_mpo, Vtr,epsilon,delta_t,n)
N= length(A);


if n < N
    if n ==1
        sleft = diag(sqrt(diag(left{n})));
        sirigth =sqrt(pinv(rigth{n}, epsilon));
        new_x =  scon({sleft,A{1},    Vtr{1}, sirigth, O{1},   left_mpo{1},rigth_mpo{1}  },...
            {[8,-1],[8,5,6],[-2,1,4] , [2, 1], [7,3,6,4], [-3,7],[5,3,2]});
        
    else
        sileft = diag(sqrt(diag(pinv(left{n}, epsilon))));
        sirigth =diag(sqrt(diag(pinv(rigth{n}, epsilon))));
        new_x = scon({A{n},    Vtr{n}, sirigth, O{n},        rigth_mpo{n}, left_mpo{n}, sileft},...
            {[7,5,6],[-2,1,4], [2,1], [8,3,6,4],[5,3,2],[7,8,9],[9,-1]});
        
        
    end  
        NA = A{n} - scon({sqrt(pinv(left{n})), new_x, conj(Vtr{n}), sqrt( pinv(rigth{n}))},...
            {[1,-1],[1,2], [2,3,-3],[-2,3]})*delta_t;
    
else
    NA =A{N};
end
end