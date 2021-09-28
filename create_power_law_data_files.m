for a=[3]
  [alpha,beta,x,fun,res] = expon_expan(a,50,6,1.e-14);
    save(['a',num2str(a),'_rcut',num2str(rcut),'_p',num2str(p), '_exp_param.mat'],'alpha','beta','X','Y');  
end