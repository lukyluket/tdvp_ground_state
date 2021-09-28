clearvars
program ='gs_tdvp';
thomas=1;
if isdeployed
    curved = str2num(getenv('curved'));
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
    d_ite_tresh =str2num(getenv('de_tresh')) ;
    d_ite=0;
    num_cent=str2num(getenv('num_cent')) ;
else
curved =0;
if curved
    chi_max=20;
    num_cent=2;
    archiv_dir  = 'effective_hamiltonian/';
    working_dir = 'saved_data/';
    name_ham =[ archiv_dir, 'eff_ham_lam_10_N_324_chi_4.mat']
    load(name_ham);
    for k=1:num_cent
    H_c{k}=Ham_eff_cent;
    end
    [chi, A,left, rigth, O, left_mpo_bg, rigth_mpo_bg, c_mpo_left, c_mpo_rigth]=...
        from_ham_2_mpo([Ham_eff_l(end:(-1):1) H_c Ham_eff_r],chi_max);
   %left{end+1}=1;
    model = 'ising';
    N=length(A);
    int_write=10;
    delta_t=0.01;
    epsilon = 1.e-14;
    delta_t = 0.01;
    t_tot = 1000;
    de_tresh =1.e-11;
    d_ite_tresh =100;
    d_ite=0;
    
    %effective_hamiltonian/eff_ham_lam_10_N_324.mat
    %c_mpo_left(1)=1;
    %ene=[1,1];
    %eneb=[];
    %de=[];
    %check = 0;
else
    if thomas
        theta=0.628;
        working_dir = 'saved_data/';
     model = 'ising_lr';
     %previous_state='gs_effective_tdvp_ising_crit_N_20_chi_10_delta_t_0.1_n_ite_1000.mat'
        init=0;
     min_ite=150;
     int_write=100;
    d=2;
    %chi_max = 40;
     

    sx = [0 ,1;1,0];
    sz =[1,0;0,-1];
    %ham =- h_xx*kron(sx,sx)-h_z/2*kron(sz,eye(2))-h_z/2*kron(eye(2),sz);
    %N = 40;
    %chi =  compute_chi_n(N,d, chi_max);
    epsilon = 1.e-14;
    delta_t = 0.01;
    t_tot = 1000;
    de_tresh =1.e-12;
    d_ite_tresh =2;
    d_ite=0;
%     %d_exp=load('/home/ltagliacozzo/Thomas/a0.1_rcut50_p9_exp_param.mat')
%     %[l_c_mpo_left{1}, l_c_mpo_rigth{1}, l_O{1}] = build_mpo(diag([1,d_exp.beta(1),1]),...
%             [0,sin(theta)*d_exp.alpha(1),0; 0,0,1;0,0,0],...
%             zeros(3),...
%             cos(theta)*[0,0,1;0,0,0;0,0,0],...
%            1,...
%             N);
%        for k_exp = 2:length(d_exp.beta)
%              [l_c_mpo_left{k_exp}, l_c_mpo_rigth{k_exp}, l_O{k_exp}] = build_mpo(diag([1,d_exp.beta(k_exp),1]),...
%             [0,sin(theta)*d_exp.alpha(k_exp),0; 0,0,1;0,0,0],...
%             zeros(3),...
%             cos(theta)*[0,0,1;0,0,0;0,0,0],...
%             1,...
%             N);
%         end
%         
%         [c_mpo_left, c_mpo_rigth, O] = concatenate_mpo(l_c_mpo_left,l_c_mpo_rigth, l_O ); 
    else
     working_dir = 'saved_data/';
     model = 'ising_lr';
     previous_state='gs_tdvp_ising_lr_N_100_chi_100_delta_t_0.01_theta_0.9891_alpha_2.mat'
     init=0;
     min_ite=150;
     int_write=100;
    d=2;
    chi_max = 200;
    sx = [0 ,1;1,0];
    sz =[1,0;0,-1];
    %ham =- h_xx*kron(sx,sx)-h_z/2*kron(sz,eye(2))-h_z/2*kron(eye(2),sz);
    N = 100;
    epsilon = 1.e-14;
    delta_t = 0.1;
    t_tot = 1000;
    de_tresh =1.e-13;
    d_ite_tresh =1.e4;
    d_ite=0;
    
    %power_law_exponent=-2;
    %num_term_approx_exp=2;
   
    %%%%%%%%%%%%INIZIALIZZO
    %left{1}=eye(chi(1));
    chi =  compute_chi_n(N,d, chi_max);
   
    %chi = chi_max* ones(1,N+1);
    %exponent = -2;
    %[prefactor, lambda]=trova_opt_exp_sum(power_law_exponent,2,num_term_approx_exp);
    %[mat_id, mat_s]=matrici_per_mpo_long_range(prefactor, lambda);
    coup_c =-1;
    coup_c(2) =-1;
    %  [l_c_mpo_left{1}, l_c_mpo_rigth{1}, l_O{1}] = build_mpo(mat_id,...
    %      mat_s,...
    %      zeros(size(mat_s)),...
    %      zeros(size(mat_s)),...
    %      1,...
    %      N);
    %
    %  [l_c_mpo_left{2}, l_c_mpo_rigth{2}, l_O{2}] = build_mpo(mat_id,...
    %      zeros(size(mat_s)),...
    %      mat_s,...
    %      zeros(size(mat_s)),...
    %      1,...
    %      N);
    % [l_c_mpo_left{3}, l_c_mpo_rigth{3}, l_O{3}] = build_mpo(mat_id,...
    %      zeros(size(mat_s)),...
    %      zeros(size(mat_s)),...
    %      mat_s,...
    %      1,...
    %      N);
    
    %   [l_c_mpo_left{1}, l_c_mpo_rigth{1}, l_O{1}] = build_mpo(eye(2),...
    %       zeros(2),...
    %       zeros(2),...
    %       [0,1;0,0],...
    %       -cos(theta),...
    %       N);
    %
    %   [l_c_mpo_left{2}, l_c_mpo_rigth{2}, l_O{2}] = build_mpo(diag([1,1]),...
     %      [0,1; 0,0],...
     %      zeros(2),...
     %      zeros(2),...
     %      -1,...
      %     N);
    
   % [l_c_mpo_left{1}, l_c_mpo_rigth{1}, l_O{1}] = build_mpo(diag([1,0,1]),...
    %    [0,1,0; 0,0,1;0,0,0],...
     %   zeros(3),...
      %  1*[0,0,1;0,0,0;0,0,0],...
       % -1,...
       % N);
    
   % [c_mpo_left, c_mpo_rigth, O] = concatenate_mpo(l_c_mpo_left,l_c_mpo_rigth, l_O );
    
    ham = construct_ising_ham(N, N+1, pi/4);
    chi_max;
    end
    %[chi, trash,trash1,trash2, O, left_mpo, rigth_mpo, c_mpo_left, c_mpo_rigth]=...
    %from_ham_2_mpo(ham,chi_max);
if init==1
    for k =1:N
        
         A{k}= rand(chi(k+1)*d, chi(k+1)*d)+1i*rand(chi(k+1)*d,chi(k+1)*d);
         [u,e] =eig( A{k});
         A{k}=reshape(u(:,1:chi(k)),chi(k+1), d,chi(k) );
          A{k}=permute(A{k},[3,1,2]);
        
%          for kk=1:d
%              A{k}(:,:,kk)=(rand(chi(k), chi(k+1))+1i*rand(chi(k),chi(k+1)))/(chi(k)*d);
%              A{k}(:,:,kk) = zeros(chi(k),chi(k+1));
%              A{k}(:,:,kk)=eye(chi(k),chi(k+1))/chi(k);
%        end
        %A{k}(1,1,1)=1;
    
       
    end
 ene = [1,1];
 eneb = [];
 de = 1;
    
    
 elseif init==0 %expand chi
     if thomas
            m_thomas = load('/home/ltagliacozzo/Thomas/data/N40_Dmax40_a0.3_theta0.628_dat.mat');
           A=m_thomas.A;
           N=length(A);
           chi_max=size(A{N/2},1);
           
           chi =  compute_chi_n(N,d, chi_max);
           theta=m_thomas.theta;
           alpha=m_thomas.a;
           %chimax=
           ene = [1,1];
           O=m_thomas.H;
           O{1}=O{2};
           O{end}=O{end-1};
           c_mpo_left=zeros(1,size(O{1},1));
           c_mpo_rigth=zeros(size(O{1},1),1);
           c_mpo_left(1)=1;
           c_mpo_left(end-1)=1;
           c_mpo_rigth(end-2)=1;
           c_mpo_rigth(end)=1;
 eneb = [];
 de = 1;
     else
     nt=0;
load([working_dir previous_state], 'A','left','rigth','ene','de','O','c_mpo_left', 'c_mpo_rigth','a' );
prev_chi=size(A{N/2},2);
left = compute_left(A,{eye(chi(1))/chi(1)});
    rigth{N} =eye(chi(N+1));

    %left=compute_left(A,left);
rigth=compute_rigth(A,rigth);
disp(' ###### EXPANDING CHI ##############')
  m_thomas = load('/home/ltagliacozzo/Thomas/data/N100_Dmax100_a2_theta0.9891_dat.mat');
   theta=m_thomas.theta;
     alpha=m_thomas.a;

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
     end

elseif init==2 %expand size
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
%     for k=1:N
%     chi_r=size(A{k},2);
%     chi_l=size(A{k},1);
%     rigth{k}=eye(chi_r);
%     left{k}=eye(chi_l)/chi_l;
%     end
  %  left = compute_left(A);
    
   % rigth{N} = eye(size(A{N},2))/size(A{N},2);
    %rigth = compute_rigth(A,rigth);
    % for k = 1: N-1
    %     H{k} = reshape( ham, d, d, d, d);
    % end
    % H{1} = H{1} -reshape (h_z*0.5*kron(eye(2),sz),d,d,d,d);
    % H{N-1} = H{N-1} -reshape( h_z*0.5*kron(sz,eye(2)),d,d,d,d);
    
    %%CHECK HAM
    %scon( {A{n}, A{n+1},  conj(A{n+1}), rigth{n+1}, H{n}, ...
    %conj(A{n}), left{n}},{[2,3,6],[3,4,5],[10,9,8],[4,9],[6,5,7,8],[11,10,7],[2,11]})
    % scon({A{N},conj(A{N}),left{N},rigth{N}},{[1,2,4],[5,3,4],[1,5],[2,3]})
    
    left = compute_left(A,{eye(chi(1))/chi(1)});
    rigth{N} =eye(chi(N+1));
    rigth = compute_rigth(A,rigth);
end
end

 
 name_file = [program '_' model '_' 'N_' num2str(N) '_chi_' ...
        num2str(chi_max) '_delta_t_' num2str(delta_t) '_theta_' num2str(theta), '_alpha_',num2str(alpha)];
%[A,left,rigth] = fix_rigth_to_identity (A, left, rigth, epsilon);
left = compute_left(A,{eye(chi(1))/chi(1)});
left_mpo = compute_left_mpo(A, O, c_mpo_left,left);
    rigth_mpo= compute_rigth_mpo(A, O, c_mpo_rigth);
    
    
    check = 0;
for nt = 1 : t_tot
    if abs(de(end)) < epsilon && nt > min_ite 
        break
    end
    d_ite=d_ite+1;
    %%COMPUTE ENERGY
    if mod(nt-1,int_write) ==0
        disp('Check');
        if fopen( [working_dir name_file '_tmp.mat']) ~= -1
            movefile( [working_dir name_file '_tmp.mat'], [working_dir 'old_' name_file '_tmp.mat'] );
        end
        save( [working_dir name_file '_tmp.mat'], 'A','O','c_mpo_left', 'c_mpo_rigth', 'ene', 'de','theta');
        
    end
    %rigth = compute_rigth(A,rigth);
   % left =  compute_left(A,{eye(chi(1))/chi(1)} );
%     if mod(nt,100) ==0
%         disp('Check');
%     end
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
    if ~curved || nt >1
     %   [A,left,rigth] = fix_rigth_to_identity (A, left, rigth, epsilon);
      %  left = compute_left(A,{eye(chi(1))/chi(1)});
       % [A, left, rigth] = fix_left_to_diag(A, left, rigth);
        %rigth = compute_rigth(A,rigth);
    end
    %  Kb = compute_kappa(A, rigth, H);
    %   norma = scon({left{N+1},rigth{N}},{[1,2],[1,2]});
    %   eneb(end+1) = trace(Kb{1}*left{1})/norma;
    %disp (['Energy_____bf g fix____'  num2str(real(trace(Kb{1}*left{1}))/norma,15) '_________Iteration_' num2str(nt)]);
    %%GO To the correct gauge
    %  v_left =left;
    %  v_A  =A;
    %  v_rigth = rigth;
    %  [A, left, rigth ] = fix_left_to_identity(A, left, rigth,epsilon);
    % rigth = compute_rigth(A,rigth);
    % [A, rigth] = fix_rigth_to_diag(A, rigth);
    % n_left = compute_left(A);
    
    %left_mpo = compute_left_mpo(A, O, c_mpo_left,left);
    rigth_mpo= compute_rigth_mpo_l(A, O, c_mpo_rigth,rigth_mpo,N-1);
    normag = scon({left{N+1},rigth{N}},{[1,2],[1,2]});
    ene(end+1) = scon({left_mpo{N},rigth_mpo{N-1}},{[1,2,3],[1,2,3]})/normag ;
    de(end+1)=ene(end-1)-ene(end);
    disp (['Energy_________'  num2str(real(ene(end)),15)....
        '____Delta_ene___',num2str(real(de(end))),'_________Iteration_', num2str(nt)]);
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
if fopen( [working_dir name_file '_tmp.mat']) ~= -1
    movefile( [working_dir name_file '_tmp.mat'], [working_dir 'old_' name_file '_tmp.mat'] );
end
save( [working_dir name_file '.mat'], 'A','O','c_mpo_left', 'c_mpo_rigth', 'ene', 'de','theta', 'alpha');
