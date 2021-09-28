clearvars
program ='gs_tdvp';
if isdeployed
    init=str2num(getenv('init'));
    chi_max=str2num(getenv('chi_max'));
    archiv_dir  = getenv( 'archiv_dir');
    working_dir = getenv('working_dir');
    ham_file =  getenv('ham_file');
    name_ham =[ archiv_dir, 'eff_ham_lam_10_N_324_chi_4.mat']
    model = getenv('model');
    epsilon = str2num(getenv('epsilon'));
    delta_t = str2num(getenv('delta_t'));
    t_tot = str2num(getenv('t_tot'));
    de_tresh = str2num(getenv('de_tresh'));
    d_ite_tresh=str2num(getenv('d_ite_tresh'));
    theta=str2num(getenv('theta')) ;
    name_long_range=getenv('name_long_range') ;
    prev_chi_max=str2num(getenv('prev_chi'));
    prev_delta_t=str2num(getenv('prev_delta_t'));
    N=str2num(getenv('N'));
    d=str2num(getenv('d'));
    min_ite=str2num(getenv('min_ite'));
    int_write=str2num(getenv('int_write'));
    alpha=str2num(getenv('alpha'));
    comp_cores=str2num(getenv('comp_cores'));
    lambda=str2num(getenv('lambda'))
    
else
    %directory where out-put is saved
     working_dir = '../saved_data/';
     mkdir(working_dir);
     %name of the model, will be part of the file-name
     model = 'ising';
     %alpha=3;
     %lambda=0.075;
     %name_long_range = 'thomas/a3_rcut50_p9_exp_param.mat';
     prev_chi_max=4;
     prev_delta_t=0.1;
     prev_theta=pi/5;
     % in case one want to use a previous state
     %previous_state='old_gs_tdvp_heis_pow_N_128_chi_128_delta_t_0.01_lambda_0.1_alpha_3_tmp.mat'
    % gs_tdvp_heis_pow_N_128_chi_128_delta_t_0.01_lambda_0.1_alpha_3_tmp.mat'
     % init = 1 start from scratch with RANDOM initial vector
     % init = 0 start from previous_state and expand chi
     % init = 2 start from previous state and expand L
     init=1;
     min_ite=150;
     int_write=10;
    d=2;
    chi_max = 4;
    sx = [0 ,1;1,0];
    sz =[1,0;0,-1];
    %ham =- h_zz*kron(sx,sx)-h_z/2*kron(sz,eye(2))-h_z/2*kron(eye(2),sz);
    N = 10;
    epsilon = 0.5e-13;
    delta_t = 0.01;
    t_tot = 100000;
    de_tresh =-1.e-13;
    d_ite_tresh =1.e4;
    d_ite = 0;
    theta = pi/5; 
    comp_cores=1;
end
maxNumCompThreads(comp_cores)

    %%%%%%%%%%%%INIZIALIZZO
    %left{1}=eye(chi(1));
    chi =  compute_chi_n(N,d, chi_max);
       d_ite=0;
%for various LR models uncomment this part anc check the files to
%understand how to build the hamiltonian as an mpo
 
%%%%%
%%%%%

% d_exp=load(name_long_range)
% if strncmp(model, 'heis' ,length(model))
%     % parameter in the Ham
%     [l_c_mpo_left, l_c_mpo_rigth, l_O] = build_hamiltonian_heis(-1,N);
%     %%xxz Hamiltonian
% elseif  strncmp(model, 'heis_pow',length(model))
%      [l_c_mpo_left, l_c_mpo_rigth, l_O] = build_hamiltonian_appr_pow_heis(d_exp,N,1);
% elseif strncmp(model, 'heis_pow_exp',length(model))
%          [l_c_mpo_left, l_c_mpo_rigth, l_O] = build_hamiltonian_appr_pow_heis(d_exp,N,lambda);
%          [l_c_mpo_left_a, l_c_mpo_rigth_a, l_O_a] = build_hamiltonian_xx_yy_zz_exp_2(lambda,N,2);
%          l_c_mpo_left{end+1}= l_c_mpo_left_a{1};
%          l_c_mpo_rigth{end+1}= l_c_mpo_rigth_a{1};
%          l_O{end+1}= l_O_a{1};
%    
% end
%         [c_mpo_left, c_mpo_rigth, O] = concatenate_mpo(l_c_mpo_left,l_c_mpo_rigth, l_O ); 
%%%%%%%
%%%%%%%
    % this is the file to change to change model
%    ham = construct_xxz_ham(N, N+1,theta);
    %ham = construct_ising_ham_p_field(N, N+1,1/4, 3/4);
  
    %%%%%
    %%%%%
    %%%Indirect version comment out if you directly use MPOS
   
    ham = construct_ising_ham(N, N+1,theta);
   
   [chi, trash,trash1,trash2, O, left_mpo, rigth_mpo, c_mpo_left, c_mpo_rigth]=...
   from_ham_2_mpo(ham,chi_max);
   %%%%%%%
    %%%%%%%


if init==1 %random start
    for k =1:N
        
         A{k}= rand(chi(k+1)*d, chi(k+1)*d)+1i*rand(chi(k+1)*d,chi(k+1)*d);
         [u,e] =eig( A{k});
         A{k}=reshape(u(:,1:chi(k)),chi(k+1), d,chi(k) );
          A{k}=permute(A{k},[3,1,2]);
             
    end
 ene = [1,1];
 eneb = [];
 de = 1;
    
    
 elseif init==0 %start from previous and expand chi
% previous_state = [program '_' model '_' 'N_' num2str(N) '_chi_' ...
       % num2str(prev_chi_max) '_delta_t_' num2str(prev_delta_t) ...
        %'_lambda_' num2str(lambda)...
        %'_alpha_' num2str(alpha)];
previous_state = [program '_' model '_' 'N_' num2str(N) '_chi_' ...
        num2str(prev_chi_max) '_delta_t_' num2str(prev_delta_t) '_theta_'...
        num2str(prev_theta)];
    
display(['Starting from' 'old_' previous_state '_tmp.mat'])   
nt=0;
load([working_dir 'old_' previous_state '_tmp.mat'], 'A','left','rigth','ene','de','O','c_mpo_left', 'c_mpo_rigth','a' );
%An =A{

prev_chi=size(A{N/2},2);
left = compute_left(A,{eye(chi(1))/chi(1)});
    rigth{N} =eye(chi(N+1));

    %left=compute_left(A,left);
rigth=compute_rigth(A,rigth);
disp(' ###### EXPANDING CHI ##############')
left = compute_left(A,{eye(chi(1))/chi(1)});
left_mpo = compute_left_mpo(A, O, c_mpo_left,left);
    rigth_mpo= compute_rigth_mpo(A, O, c_mpo_rigth);
    d_chi=chi_max-prev_chi;
        NA=expand_chi_2(A,O,d_chi,left,rigth,left_mpo, rigth_mpo,c_mpo_rigth, epsilon,delta_t);
        A=NA;
        rigth = compute_rigth(A,rigth);
        left =  compute_left(A, left);
        [A,left,rigth] = fix_rigth_to_identity (A, left, rigth, epsilon);
        left = compute_left(A,{eye(chi(1))/chi(1)});
        [A, left, rigth] = fix_left_to_diag(A, left, rigth);
        left_mpo = compute_left_mpo(A, O, c_mpo_left,left);
        rigth_mpo= compute_rigth_mpo(A, O, c_mpo_rigth);
        normag = scon({left{N+1},rigth{N}},{[1,2],[1,2]});
        ene (end+1)= scon({left_mpo{2},rigth_mpo{1}},{[1,2,3],[1,2,3]})/normag ;
        de(end+1)=ene(end-1)-ene(end);
        disp (['Energy_________'  num2str(real(ene(end)),15)....
            '____Delta_ene___',num2str(real(de(end))),'_________Iteration_', num2str(nt)]);

elseif init==2 % start from previous and expand N
  load([working_dir previous_state], 'A','left','rigth','ene','de');
  disp(' ###### EXPANDING N ##############')
%left = compute_left(A,{eye(chi(1))/chi(1)});

  previous_N=length(A);
  diff_N= N-previous_N;
  pos=previous_N/2*ones(1,diff_N);
  NA(1:previous_N/2) = A(1: previous_N/2);
  NA(previous_N/2+1:previous_N/2+diff_N)=A(pos);
  NA(previous_N/2+diff_N+1:previous_N/2+diff_N+previous_N/2)=A(previous_N/2+1:end);
  A=NA;
  rigth{N} =eye(chi(N+1));
  rigth=compute_rigth(A,rigth);
  left=compute_left(A,{eye(chi(1))/chi(1)});
     
end

% Now everything is setteled for security we recompute left and rigth


left = compute_left(A,{eye(chi(1))/chi(1)});
rigth{N} =eye(chi(N+1));
rigth = compute_rigth(A,rigth);

 
 %name_file = [program '_' model '_' 'N_' num2str(N) '_chi_' ...
  %      num2str(chi_max) '_delta_t_' num2str(delta_t) '_lambda_'...
   %     num2str(lambda) '_alpha_' num2str(alpha)];
 name_file = [program '_' model '_' 'N_' num2str(N) '_chi_' ...
        num2str(chi_max) '_delta_t_' num2str(delta_t) '_theta_'...
        num2str(theta)];

%[A,left,rigth] = fix_rigth_to_identity (A, left, rigth, epsilon);
%left = compute_left(A,{eye(chi(1))/chi(1)});
left_mpo = compute_left_mpo(A, O, c_mpo_left,left);
rigth_mpo= compute_rigth_mpo(A, O, c_mpo_rigth);
    
    
    check = 0;
%here we start optimization, loop at most up to t_tot, break if de< epsilon
%and at least min_ite have been done (important if one starts from a ground
%state similar to the one desired, optimization can be slow at the
%beginning.
 for nt = 1 : t_tot
    if abs(de(end)) < epsilon && nt > min_ite 
        break
    end
    d_ite=d_ite+1;
    %%COMPUTE ENERGY
    %write to file every int_write iterations in case the simulation get
    %killed one can restart from the files _tmp.mat
    if mod(nt-1,int_write) ==0
        disp('Check');
        d_f= fopen( [working_dir name_file '_tmp.mat']);
        if  d_f~= -1
            fclose(d_f);
            movefile( [working_dir name_file '_tmp.mat'], [working_dir 'old_' name_file '_tmp.mat'] );
        end
        save( [working_dir name_file '_tmp.mat'], 'A','O','c_mpo_left', 'c_mpo_rigth', 'ene', 'de','theta');
        
    end
    if check
        left_mpo_bg = compute_left_mpo(A, O, c_mpo_left,left);
        rigth_mpo_bg = compute_rigth_mpo(A, O, c_mpo_rigth);
        norma = scon({left{N+1},rigth{N}},{[1,2],[1,2]});
        eneb(end+1) = scon({left_mpo_bg{3},rigth_mpo_bg{2}},{[1,2,3],[1,2,3]})/norma;
        v_left =left;
        v_rigth = rigth;
        v_A =A;
        disp (['Energy_____bf g fix____'  num2str(real(eneb(end)),15) '_________Iteration_' num2str(nt)]);
    end
   
    rigth_mpo= compute_rigth_mpo_l(A, O, c_mpo_rigth,rigth_mpo,N-1);
    normag = scon({left{N+1},rigth{N}},{[1,2],[1,2]});
    ene(end+1) = scon({left_mpo{N},rigth_mpo{N-1}},{[1,2,3],[1,2,3]})/normag ;
    de(end+1)=ene(end-1)-ene(end);
    disp (['Energy_________'  num2str(real(ene(end)),15)....
        '____Delta_ene___',num2str(real(de(end))),'_________Iteration_', num2str(nt)]);
   %one can decide to automatically expand chi, however i tend to avoid
   %this by putting de_tresh to a neg. value
    if de(end) < de_tresh && d_ite > d_ite_tresh
        d_ite=0;
        disp(' ###### EXPANDING CHI ##############')
        NA=expand_chi_2(A,O,chi_max,left,rigth,left_mpo, rigth_mpo,c_mpo_rigth, epsilon,delta_t);
        A=NA;
        rigth = compute_rigth(A,rigth);
        left =  compute_left(A, left);
        chi_max=size(A{N/2},1);
        chi=compute_chi_n(N, d, chi_max);
         name_file = [program '_' model '_' 'N_' num2str(N) '_chi_' ...
        num2str(chi_max) '_delta_t_' num2str(delta_t) '_theta_' num2str(theta), '_alpha_',num2str(alpha)];

        [A,left,rigth] = fix_rigth_to_identity (A, left, rigth, epsilon);
        left = compute_left(A,{eye(chi(1))/chi(1)});
        [A, left, rigth] = fix_left_to_diag(A, left, rigth);
        left_mpo = compute_left_mpo(A, O, c_mpo_left,left);
        rigth_mpo= compute_rigth_mpo(A, O, c_mpo_rigth);
        normag = scon({left{N+1},rigth{N}},{[1,2],[1,2]});
        ene(end+1) = scon({left_mpo{2},rigth_mpo{1}},{[1,2,3],[1,2,3]})/normag ;
        de(end+1)=ene(end-1)-ene(end);
        disp (['Energy_________'  num2str(real(ene(end)),15)....
            '____Delta_ene___',num2str(real(de(end))),'_________Iteration_', num2str(nt)]);
    end
    %start from  rigth sweep of optimization
    for n = N :(-1) :1
        
        Vtr = fix_gauge_for_B_l (A, rigth,n);
        A{n} = update_A_mpo_l(A, O,left, rigth, left_mpo, rigth_mpo, Vtr,epsilon,delta_t,n);
        %rigth = compute_rigth_l(A,rigth,n-1);
        [A,left,rigth] = fix_rigth_to_identity_l (A, left, rigth, epsilon,n-1);
        if n > 1
            rigth_mpo= compute_rigth_mpo_l(A, O, c_mpo_rigth,rigth_mpo,n-1);
        end
        
        
        %left = compute_left_l(A,left);
    end
    [A,left,rigth] = fix_left_to_diag_l (A, {eye(chi(1))/chi(1)} , rigth, 1);
        left_mpo= compute_left_mpo_l( A, O, left_mpo, c_mpo_left,left,1);
    %start from left sweep of optimization
        for n = 1 :N
       
        
       % left_mpo= compute_left_mpo_l( A, O, left_mpo, cl_O_left,n);
        
        
        Vtr = fix_gauge_for_B_l (A, rigth,n);
        A{n} = update_A_mpo_l(A, O,left, rigth, left_mpo, rigth_mpo, Vtr,epsilon,delta_t,n);
        %rigth = compute_rigth_l(A,rigth,n-1);
        [A,left,rigth] = fix_left_to_diag_l (A, left, rigth, n+1);
        if n < N
            left_mpo= compute_left_mpo_l(A, O,left_mpo,c_mpo_left,left ,n+1);
        end
        
        
        %left = compute_left_l(A,left);
    end
    
    %A =NA;
end
%energies = compute_energy(A,left,rigth,H);
d_f = fopen( [working_dir name_file '_tmp.mat']);
if  d_f~= -1
    fclose(d_f);
    movefile( [working_dir name_file '_tmp.mat'], [working_dir 'old_' name_file '_tmp.mat'] );
end
%save the data at the end
save( [working_dir name_file '.mat'], 'A','O','c_mpo_left', 'c_mpo_rigth', 'ene', 'de','theta');
