function [alpha,beta,x,func,res] = expon_expan(alph,rcut,N,tol,previous)

clear vars
prev_length=length(previous);
exitfl=0;
res=1;
beta = 10^(rcut/N);
while exitfl == 0 || res >= tol   %e-10
    x0=rand(2*N,1);
    x0(1:prev_length) = previous;
%     for i = 0: rcut-1
%    lambda(i+1) =  -alph /(beta ^i);
%    prefactor(i+1) = exp(-alph)*beta ^(-i*alph);
% end
% 
% overall_constant = exp(lambda)*prefactor.';
% prefactor = prefactor./overall_constant;

    
    options = optimset('MaxFunEvals',N*10000,'MaxIter',10000,'TolFun',tol);
    [coef,res,~,exitfl] = lsqnonlin(@lsq_exp,x0);
    res
end
disp(res);
alpha=coef(1:2:end);
beta=coef(2:2:end);

for r= 1:rcut
    func(r)=0;
    for i= 1:N
        func(r)=func(r)+alpha(i)*beta(i)^(r-1);
    end
end

x=1:rcut;

function [y] = lsq_exp(x)
    y=zeros(rcut);
    for r = 1:rcut
        for i=1:N
            y(r)=y(r)+x(2*i-1)*x(2*i)^(r-1);
        end
        y(r) =y(r) -1/r^alph;
    end
end
figure
loglog(x,func,'.')
hold all
loglog(x ,1./x.^alph, 'O')
end
