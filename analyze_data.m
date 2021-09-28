function analyze_data()
mkdir('num_plots')
repository = '/home/luca/data/saved_data/'
ent_fig=figure;
fig_corr=figure;
%fig_ent_h =figure;
% ren_ent_fig_2 =figure;
% ren_ent_fig_3 =figure;
% ren_ent_fig_4 =figure;
% ren_ent_fig_5 =figure;
lambda=0.1
n_mod=2;

%for chi = [32 64 128]
for n_mod=[ 2 ]
    if n_mod ==0
        model='heis';
        model_l = 'Heis';
    elseif n_mod==1
        model='heis_pow';
        model_l ='Heis-Power'
    elseif n_mod==2
        model='heis_pow_exp';
        model_l ='Heis-Power-exp';
    end   
ent_half_chain=[];
ent_half_chain_p=[];
ent_half_chain_m=[];
lunghezza=[];
corr_xx=[];
corr_yy=[];
corr_zz=[];
s_x=[0,1;1,0]
s_z=diag([1,-1]);
%N='8'
alpha='3'
epsilon =1.e-13;

%corr_fig=figure;
%array_files=ls([repository 'gs_tdvp_' model '_N*' 'chi_' num2str(chi) '*' 'alpha_' alpha '_tmp.mat']);
%list_files =strsplit(array_files);
if n_mod==0
%list_files{1}=[ repository 'gs_tdvp_heis_N_8_chi_64_delta_t_0.01_lambda_0.1_alpha_3.mat'];
%list_files{2}=[ repository 'gs_tdvp_heis_N_16_chi_64_delta_t_0.01_lambda_0.1_alpha_3.mat'];
%list_files{3}=[ repository 'gs_tdvp_heis_N_32_chi_64_delta_t_0.01_lambda_0.1_alpha_3.mat'];
%list_files{1}=[ repository 'gs_tdvp_heis_N_64_chi_128_delta_t_0.01_lambda_0.1_alpha_3.mat'];
list_files{1}=[ repository 'gs_tdvp_heis_N_128_chi_128_delta_t_0.01_lambda_0.1_alpha_3_tmp.mat']; %questo e' ancora temp
elseif n_mod==1
%list_files{1}=[ repository 'gs_tdvp_heis_pow_N_8_chi_64_delta_t_0.01_lambda_0.1_alpha_3.mat'];
%list_files{2}=[ repository 'gs_tdvp_heis_pow_N_16_chi_128_delta_t_0.01_lambda_0.1_alpha_3.mat'];
%list_files{3}=[ repository 'gs_tdvp_heis_pow_N_32_chi_128_delta_t_0.01_lambda_0.1_alpha_3.mat'];
%list_files{4}=[ repository 'gs_tdvp_heis_pow_N_64_chi_128_delta_t_0.01_lambda_0.1_alpha_3.mat'];
list_files{1}=[ repository 'gs_tdvp_heis_pow_N_64_chi_128_delta_t_0.01_lambda_0.1_alpha_3_tmp.mat'];
elseif n_mod==2
 %list_files{1}=[ repository 'gs_tdvp_heis_pow_exp_N_8_chi_64_delta_t_0.01_lambda_' num2str(lambda) '_alpha_3.mat'];
 %list_files{2}=[ repository 'gs_tdvp_heis_pow_exp_N_16_chi_128_delta_t_0.1_lambda_' num2str(lambda) '_alpha_3.mat'];
 %list_files{3}= [repository 'gs_tdvp_heis_pow_exp_N_32_chi_128_delta_t_0.1_lambda_' num2str(lambda)  '_alpha_3.mat'];
 %list_files{4}= [repository 'gs_tdvp_heis_pow_exp_N_64_chi_128_delta_t_0.1_lambda_0.1_alpha_3.mat'];
 %list_files{4}= [repository 'gs_tdvp_heis_pow_exp_N_64_chi_256_delta_t_0.1_lambda_' num2str(lambda) '_alpha_3_tmp.mat'];
% list_files{5}= [repository 'gs_tdvp_heis_pow_exp_N_128_chi_64_delta_t_0.1_lambda_0.1_alpha_3.mat'];
 %list_files{6}= [repository 'gs_tdvp_heis_pow_exp_N_128_chi_128_delta_t_0.1_lambda_0.1_alpha_3_tmp.mat'];
 list_files{1}= [repository 'gs_tdvp_heis_pow_exp_N_128_chi_128_delta_t_0.5_lambda_0.1_alpha_3.mat'];
 list_files{3}= [repository 'gs_tdvp_heis_pow_exp_N_128_chi_128_delta_t_0.5_lambda_0.05_alpha_3.mat'];
 list_files{2}= [repository 'gs_tdvp_heis_pow_exp_N_128_chi_128_delta_t_0.5_lambda_0.075_alpha_3.mat'];

end    
alpha=str2num(alpha);
%N=str2num(N);
list_lambdas = [0.1, 0.075, 0.05]
for n=1:length(list_files)
load(list_files{n})
N=length(A);
chi=size(A{N/2+1},2)
point_1 = N/2-N/4;
point_2 =N/2 +N/4;
lambda =list_lambdas(n);
left = compute_left(A,{1});
rigth{N} =eye(1);
rigth = compute_rigth(A,rigth);
[A,left,rigth] = fix_rigth_to_identity (A, left, rigth, epsilon);
 left = compute_left(A,{1});
[A, left, rigth] = fix_left_to_diag(A, left, rigth);
[entropies,norm_state] = compute_entropies(left);
[r2_entropies,norm_state] = compute_ren_entropies(left,2);
[r3_entropies,norm_state] = compute_ren_entropies(left,3);
[r4_entropies,norm_state] = compute_ren_entropies(left,4);
[r5_entropies,norm_state] = compute_ren_entropies(left,5);

%[dist,correlator_xx] = compute_corr(A,s_x);

ent_half_chain(end+1)=entropies(N/2);
ent_half_chain_m(end+1)=entropies(N/2-1);
ent_half_chain_p(end+1)=entropies(N/2+1);

lunghezza(end+1)=N;
%corr_xx(end+1)=two_point_corr(A,s_x, point_1,point_2)
%corr_zz(end+1)=two_point_corr(A,s_z, point_1,point_2)

%l_theta(n)=theta;

figure(ent_fig)
%h=plot(log(2*N./pi*sin(pi*(1:N-1)/N)), entropies, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
h=plot(log(2*N./pi*sin(pi*(2:2:N-1)/N)), entropies(2:2:end), 'O', 'DisplayName', [model_l ' \lambda = ' num2str(lambda) ' N = ' num2str(N)])
set(h, 'MarkerFaceColor', get(h, 'Color'));
xl=xlabel('$\log(\frac{2 \pi}{N} sin(\frac{\pi l}{N}))$')
set(xl, 'Interpreter', 'latex')
ylabel('S')
hold all
legend OFF
legend SHOW
%end
% figure(norm_fig)
% plot(real(norm_state)-1, '-O')
% 
% figure(fig_corr)
% %z11b= sin(pi*(N-dist)/N);
% %z21b= sin(pi*N/N);
% %z12 =sin(pi*dist/N);
% pos1 = N/2 -dist./2;
% pos2 = pos1 +dist;
% 
% x  = -(2-2*cos(pi*dist/N))./ 4./ sin(pi*pos1/N)./sin(pi*pos2/N);
% 
% norm_corr= (( sin(pi*pos1/N).*sin(pi*pos2/N) ))  .*correlator_xx;
% %eta= z11b.^2./z21b./z21b;
% 
% %h=plot((1./(2*cosh(pi*dist/N)-2).^0.5)*pi/N, correlator_xx, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% h=plot(x,norm_corr, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% 
% set(h, 'MarkerFaceColor', get(h, 'Color'));
% xl=xlabel('$\frac{\pi}{N} /\sqrt{cosh( \frac{\pi l}{N})-1}$')
% set(xl, 'Interpreter', 'latex')
% ylabel('C_xx')
% hold all
% legend OFF
% legend SHOW
% 
% hold on
n_start=N/2-8;
n_end=N/2 -2 ;
coeff=polyfit (log(2*N./pi*sin(pi*(n_start:2:n_end)/N)), entropies(n_start:2:n_end),1)
central_charge(n)=coeff(1)*6;

%plot(log(2*N./pi*sin(pi*(n_start:2:n_end)/N)), coeff(1)*log(2*N./pi*sin(pi*(n_start:2:n_end)/N))+coeff(2),'-k')
% figure(ren_ent_fig_2)
% %h=plot((2*N./pi*sin(pi*(1:N-1)/N)).^(1/12*(2-1/2)), r2_entropies, '-O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% h=plot(1:N-1, r2_entropies, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% 
% set(h, 'MarkerFaceColor', get(h, 'Color'));
% %xl=xlabel('$(\frac{2 \pi}{N} sin(\frac{\pi l}{N}))^{c/12(n-1/n)}$')
% xlabel('l')
% %set(xl, 'Interpreter', 'latex')
% ylabel('R_2')
% hold all
% 
% figure(ren_ent_fig_3)
% %h=plot((2*N./pi*sin(pi*(1:N-1)/N)).^(1/12*(2-1/2)), r2_entropies, '-O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% h=plot(1:N-1, r3_entropies, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% 
% set(h, 'MarkerFaceColor', get(h, 'Color'));
% %xl=xlabel('$(\frac{2 \pi}{N} sin(\frac{\pi l}{N}))^{c/12(n-1/n)}$')
% %set(xl, 'Interpreter', 'latex')
% xlabel('l')
% ylabel('R_3')
% hold all
% figure(ren_ent_fig_4)
% %h=plot((2*N./pi*sin(pi*(1:N-1)/N)).^(1/12*(2-1/2)), r2_entropies, '-O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% h=plot(1:N-1, r4_entropies, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% 
% set(h, 'MarkerFaceColor', get(h, 'Color'));
% %xl=xlabel('$(\frac{2 \pi}{N} sin(\frac{\pi l}{N}))^{c/12(n-1/n)}$')
% %set(xl, 'Interpreter', 'latex')
% xlabel('l')
% ylabel('R_4')
% hold all
% figure(ren_ent_fig_5)
% %h=plot((2*N./pi*sin(pi*(1:N-1)/N)).^(1/12*(2-1/2)), r2_entropies, '-O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% h=plot(1:N-1, r5_entropies, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
% 
% set(h, 'MarkerFaceColor', get(h, 'Color'));
% %xl=xlabel('$(\frac{2 \pi}{N} sin(\frac{\pi l}{N}))^{c/12(n-1/n)}$')
% %set(xl, 'Interpreter', 'latex')
% xlabel('l')
% ylabel('R_5')
% hold all
end
figure
pc =plot(list_lambdas, central_charge, 'O')
set(pc, 'MarkerFaceColor', get(pc, 'Color'));
x_points = 1:N-1;
%[coef,res,exitfl,f_values]= fss_sca_osc(2:2:N-1,r2_entropies(2:2:end),[0.73],N);
% [coef_2,res_2,exitfl,f_values_2]= fss_sca_osc(x_points,r2_entropies(2:2:end),[0.73],N,2,1);
% [coef_2_2,res_2,exitfl,f_values_21]= fss_sca_osc(x_points,r2_entropies(2:2:end)./f_values_2,[0.01,-1/4],N,2,0,coef_2);
% 
% [coef_3,res_3,exitfl,f_values_3]= fss_sca_osc(x_points,r3_entropies(2:2:end),[0.73],N,3,1);
% [coef_3_2,res_3,exitfl,f_values_31]= fss_sca_osc(x_points,r3_entropies(2:2:end)./f_values_3,[0.01,-1/4],N,3,0,coef_3);
% 
%[coef,res,exitfl,f_values]= fss_sca_osc(2:2:N-1,r2_entropies,[0.73],N);

load /home/luca/dismissed_dropbox/3blo/ren_2_heis_chi_8_N_30.mat
x_points=1:max_block_size;
N=30
figure(ren_ent_fig_2)
hold off
h =plot(1:max_block_size,reny1_2, 'o')
%set(h, 'MarkerFaceColor', get(h, 'Color'));

[coef_2,res_2,exitfl,f_values_2]= fss_sca_osc(x_points,reny1_2,[0.73],N,2,1);
 [coef_2_2,res_2,exitfl,f_values_21]= fss_sca_osc(x_points,reny1_2./f_values_2,[1,-1/4,1],N,2,0,coef_2);


 %[coef_2,res_2,exitfl,f_values_2]= fss_sca_osc(x_points,r2_entropies,[0.73],N,2,1);
 %[coef_2_2,res_2,exitfl,f_values_21]= fss_sca_osc(x_points,r2_entropies./f_values_2,[0.01,-1/4,1],N,2,0,coef_2);
 
%  [coef_3,res_3,exitfl,f_values_3]= fss_sca_osc(x_points,r3_entropies,[0.73],N,3,1);
%  [coef_3_2,res_3,exitfl,f_values_31]= fss_sca_osc(x_points,r3_entropies./f_values_3,[0.01,-1/6,1],N,3,0,coef_3);
%  
%  
%  [coef_4,res_4,exitfl,f_values_4]= fss_sca_osc(x_points,r4_entropies,[0.73],N,4,1);
%  [coef_4_2,res_4,exitfl,f_values_41]= fss_sca_osc(x_points,r4_entropies./f_values_4,[0.01,-1/8,1],N,4,0,coef_4);
% % 
% 
% [coef_5,res_5,exitfl,f_values_5]= fss_sca_osc(x_points,r5_entropies,[0.73],N,5,1);
%  [coef_5_2,res_5,exitfl,f_values_51]= fss_sca_osc(x_points,r5_entropies./f_values_5,[0.01,-1/10,1],N,5,0,coef_5);

% [coef_3,res_3,exitfl,f_values_3]= fss_sca_osc(x_points,r3_entropies(2:2:end),[0.73,0.01,-1/6],N,3);
% [coef_4,res_2,exitfl,f_values_4]= fss_sca_osc(x_points,r4_entropies(2:2:end),[0.73,0.01,-1/8],N,4);
% [coef_5,res_3,exitfl,f_values_5]= fss_sca_osc(x_points,r5_entropies(2:2:end),[0.73,0.01,-1/10],N,5);

figure(ren_ent_fig_2)
hold all
%plot(2:2:N-1,f_values, '.')
plot(x_points,f_values_2, '-<')


figure
%h=plot(x_points, r2_entropies(2:2:end)./f_values_2, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
%h=plot(x_points, r2_entropies./f_values_2, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
h=plot(x_points, reny1_2./f_values_2, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])

set(h, 'MarkerFaceColor', get(h, 'Color'));
hold all
plot(x_points,f_values_21,'-<')
xlabel('l')
ylabel('R2/CFT')
figure(ren_ent_fig_3)
plot(x_points,f_values_3, '-<')


figure 

%h=plot(x_points, r3_entropies(2:2:end)./f_values_3, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
h=plot(x_points, r3_entropies./f_values_3, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])

set(h, 'MarkerFaceColor', get(h, 'Color'));
hold all
plot(x_points,f_values_31,'-<')
xlabel('l')
ylabel('R3/CFT')


figure(ren_ent_fig_4)

%plot(2:2:N-1,f_values, '.')
plot(x_points,f_values_4, '-<')

figure 

%h=plot(x_points, r3_entropies(2:2:end)./f_values_3, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
h=plot(x_points, r4_entropies./f_values_4, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])

set(h, 'MarkerFaceColor', get(h, 'Color'));
hold all
plot(x_points,f_values_41,'-<')
xlabel('l')
ylabel('R4/CFT')




figure(ren_ent_fig_5)
plot(x_points,f_values_5, '-<')


figure 

%h=plot(x_points, r3_entropies(2:2:end)./f_values_3, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])
h=plot(x_points, r5_entropies./f_values_5, 'O', 'DisplayName', [model_l ' \chi = ' num2str(chi) ' N = ' num2str(N)])

set(h, 'MarkerFaceColor', get(h, 'Color'));
hold all
plot(x_points,f_values_51,'-<')
xlabel('l')
ylabel('R5/CFT')


figure(ent_fig)
nN=128;
  plot(log(2*nN/pi*sin(pi*(1:nN-1)/nN)),1/6*log(2*nN/pi*sin(pi*(1:nN-1)/nN))+0.23, 'k', 'DisplayName', 'CFT')
hold all


figure(fig_corr)
plot((1./(2*cosh(pi*dist/N)-2).^0.5)*pi/N, (log((1./(2*cosh(pi*dist/N)-2).^0.5)*pi/N+1)).^0.7, '-k', 'DisplayName', 'CFT')
[lunghezza, i_sort]= sort(lunghezza, 'ascend');
corr_xx= corr_xx(i_sort);
corr_zz= corr_zz(i_sort);
ent_half_chain = ent_half_chain(i_sort);
ent_half_chain_m =ent_half_chain_m(i_sort);



%h=plot(log(lunghezza/2), log(corr_xx), 'O-','DisplayName', ['C_{xx}-' model_l ' \chi =' num2str(chi) ])
%set(h, 'MarkerFaceColor', get(h, 'Color'));
%xlabel('log(N/2)')
%ylabel('log(corr)')
%hold all
%h=plot(log(lunghezza/2), log(corr_zz), '<-','DisplayName', ['C_{zz}-' model_l ' \chi =' num2str(chi)])
%set(h, 'MarkerFaceColor', get(h, 'Color'));
%xlabel('log(N/2)')
%ylabel('log(C)')
%legend OFF
%legend SHOW
figure(fig_ent_h)
h=plot(log(lunghezza/2), ent_half_chain, 'O-','DisplayName', ['S_{N/2}-' model_l ' \chi =' num2str(chi)])
set(h, 'MarkerFaceColor', get(h, 'Color'));

xlabel('log(N/2)')
ylabel('S')
hold all
h=plot(log(lunghezza/2), ent_half_chain_m, 'd-','DisplayName', ['S_{N/2-1}' model_l ' \chi =' num2str(chi)])
set(h, 'MarkerFaceColor', get(h, 'Color'));
legend OFF
legend SHOW
end
%end
saveas(fig_corr,['num_plots/corr_' model '_chi_' num2str(chi) '.pdf'], 'pdf')
saveas(fig_ent_h,['num_plots/ent_h_' model '_chi_' num2str(chi) '.pdf'], 'pdf')
saveas(ent_fig,['num_plots/ent_' model '_chi_'  num2str(chi) '.pdf'], 'pdf')

saveas(fig_corr,['num_plots/corr_' model '_chi_' num2str(chi) '.fig'], 'fig')
saveas(fig_ent_h,['num_plots/ent_h' model '_chi_' num2str(chi) '.fig'], 'fig')
saveas(ent_fig,['num_plots/ent_' model '_chi_' num2str(chi) '.fig'], 'fig')

saveas(fig_corr,['num_plots/corr_' model '_chi_' num2str(chi) '.jpg'], 'jpg')
saveas(fig_ent_h,['num_plots/ent_h_' model '_chi_' num2str(chi) '.jpg'], 'jpg')
saveas(ent_fig,['num_plots/ent__' model '_chi_' num2str(chi) '.jpg'], 'jpg')


return


figure
plot(l_theta, ent_half_chain, 'O', 'DisplayName', 'N/2')
hold on
plot(l_theta, ent_half_chain_m, 'O','DisplayName', 'N/2-1')
plot(l_theta, ent_half_chain_p, 'O','DisplayName', 'N/2+1')
xlabel('\theta')
ylabel('S_{N/2}')
figure
plot(l_theta, cos(l_theta)./sin(l_theta), 'O', 'MarkerFaceColor', 'b')
hold all
plot(l_theta, ones(length(l_theta)),'-r')
ylabel('\Delta')
xlabel('\theta')

end

function [entropies,norm_state]=compute_entropies(left)
for l=2:length(left)-1
schmidt_values=diag(left{l});
reg_schimdt_values = schmidt_values(find(schmidt_values >1e-13));
entropies(l-1)=-real(reg_schimdt_values.'*log(real(reg_schimdt_values)));
norm_state(l-1)=sum(reg_schimdt_values);

end
end
function [entropies,norm_state]=compute_ren_entropies(left,n)
for l=2:length(left)-1
schmidt_values=diag(left{l});
reg_schimdt_values = schmidt_values(find(schmidt_values >1e-13));
entropies(l-1)=(sum(reg_schimdt_values.^n));
norm_state(l-1)=sum(reg_schimdt_values);

end
end
function [dist,correlator]=compute_corr(A,loc_op)
correlator=[];
dist=[];
N = length(A);
for d =1:N/2-1
    point_1=N/2-d;
    point_2=N/2+d;
    dist(end+1)=2*d;
    correlator(end+1) =two_point_corr(A,loc_op, point_1,point_2);
end
end
function corr_xx=two_point_corr(A,loc_op, point_1,point_2)
N=length(A);
A_x =A;

A_x{point_1} =scon({A_x{point_1},loc_op},{[-1,-2,1],[1,-3]});
A_x{point_2} =scon({A_x{point_2},loc_op},{[-1,-2,1],[1,-3]});

left_x = compute_mix_left(A_x,A,{1});
rigth_x{N} =eye(1);
rigth_x = compute_mix_rigth(A_x,A,rigth_x);
corr_xx = scon({left_x{N+1},rigth_x{N}},{[1,2],[1,2]});
end
