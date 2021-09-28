
function [chi_n,A,left, rigth, O, left_mpo_bg, rigth_mpo_bg, ...
    c_mpo_left, c_mpo_rigth] = from_ham_2_mpo(ham,chi)
%clearvars
%L = 30;
cut =1.e-11;
%d= 2;
%chi=10;
L= length(ham);

for k=1:(2):L
    d(k) = size(ham{k},1);
    d(k+1) = size(ham{k},2);
    
end
d(L+1) = size(ham{L},2);
chi_n =  compute_chi_n(L+1,d, chi);
for k=1:(2):L
    
    for i=1:d(k)
        A{k}(:,:,i)=(rand(chi_n(k),chi_n(k+1))+1i*rand(chi_n(k),chi_n(k+1)))/(chi_n(k+1)*d(k));
        %A{k+1}(:,:,i)=rand(chi)+1i*rand(chi);
    end
    for i=1:d(k+1)
        %A{k}(:,:,i)=rand(chi)+1i*rand(chi);
        A{k+1}(:,:,i)=(rand(chi_n(k+1),chi_n(k+2))+1i*rand(chi_n(k+1),chi_n(k+2)))/(chi_n(k+1)*d(k+1));
    end
    
end
for i=1:d(L+1)
    %A{k}(:,:,i)=rand(chi)+1i*rand(chi);
    A{L+1}(:,:,i)=(rand(chi_n(L+1),chi_n(L+2))+1i*rand(chi_n(L+1),chi_n(L+2)))/(chi_n(L+1)*d(L+1));
end
%for k =1:L-1
%   ham{k} =rand(d,d,d,d)+1i*rand(d,d,d,d);
%end
rigth{L+1}=eye(chi_n(L+2));
%rigth = compute_rigth(A,rigth);
%left =  compute_left(A, {eye(chi_n(1))});
for k=1:L+1
    chi_r=size(A{k},2);
    chi_l=size(A{k},1);
    rigth{k}=eye(chi_r);
    left{k}=eye(chi_l)/chi_l;
end

Kb = compute_kappa(A, rigth, ham);
norma = scon({left{L+1},rigth{L}},{[1,2],[1,2]});
ene_ham = trace(Kb{1}*left{1})/norma;
figure

for k = 1: length(ham)
    ham_m{k} = reshape(permute(ham{k},[1,3,2,4]),d(k)*d(k),d(k+1)*d(k+1));
    [u{k},s{k},v{k}] = svd(ham_m{k});
    [uc{k},sc{k},vc{k}] = svd(ham_m{k}.');
    i_s{k}= find(diag(s{k})>cut);
    tot_size(k) =length(i_s{k});
    semilogy(diag(s{k}));
    hold all
end

dim_mpo= max(tot_size) +2;
for k =1:length(ham)
    for i_op=1:tot_size(k)
        vd= v{k}';
        op_left{i_op}=reshape(u{k}(:,i_op)*sqrt(s{k}(i_op,i_op)),d(k),d(k));
        op_rigth{i_op}=reshape(vd(i_op,:)*sqrt(s{k}(i_op,i_op)),d(k+1),d(k+1));
    end
    
    %     rec_ham =0;
    %     for i_op=1:length(i_s{k})
    %     rec_ham = rec_ham +reshape(op_left{i_op},d^2,1)*reshape(op_rigth{i_op},1,d^2);
    %     end
    %
    %     max(ham_m{k}(1:end)-rec_ham(1:end))
    O{k}(1,1,:,:) =eye(d(k));
    O{k}(dim_mpo,dim_mpo,:,:) =eye(d(k));
    for i_op=1:length(i_s{k})
        O{k}(1,i_op+1,:,:) = op_left{i_op};
        O{k+1}(i_op+1,dim_mpo,:,:) = op_rigth{i_op};
    end
    
end
O{k+1}(1,1,:,:) =eye(d(k+1));
O{k+1}(dim_mpo,dim_mpo,:,:) =eye(d(k+1));
c_mpo_left=zeros(1,dim_mpo);
c_mpo_rigth=zeros(dim_mpo,1);
c_mpo_left(1)=1;
c_mpo_rigth(end)=1;

left_mpo_bg = compute_left_mpo(A, O, c_mpo_left,left);
rigth_mpo_bg = compute_rigth_mpo(A, O, c_mpo_rigth);
norma = scon({left{L+1},rigth{L}},{[1,2],[1,2]});
ene_mpo = scon({left_mpo_bg{2},rigth_mpo_bg{1}},{[1,2,3],[1,2,3]})/norma;
ene_ham-ene_mpo
end