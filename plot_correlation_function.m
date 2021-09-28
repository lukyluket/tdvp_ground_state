function [dist,correlator]=plot_correlation_function(model, N,chi_max, delta_t, theta,op)
s_x =[0,1/2;1/2,0];
s_z=diag([1/2,-1/2]);
program='gs_tdvp';
working_dir='../saved_data/'
epsilon =1.e-13;

%N=10
%chi_max=10
%delta_t=0.1
%theta=pi/5
state = [program '_' model '_' 'N_' num2str(N) '_chi_' ...
        num2str(chi_max) '_delta_t_' num2str(delta_t) '_theta_'...
        num2str(theta)];
 load([working_dir  state '.mat'], 'A','left','rigth','ene','de','O','c_mpo_left', 'c_mpo_rigth','a' );
left = compute_left(A,{1});
rigth{N} =eye(1);
rigth = compute_rigth(A,rigth);
[A,left,rigth] = fix_rigth_to_identity (A, left, rigth, epsilon);
 left = compute_left(A,{1});
[A, left, rigth] = fix_left_to_diag(A, left, rigth);

corr_xx=two_point_corr(A,loc_op, point_1,point_2)
if op==1
 [dist,correlator]=compute_average(A,s_x)
    %[dist,correlator]=compute_corr(A,s_x)
elseif op==2
     [dist,correlator]=compute_average(A,s_z)

  [dist,correlator]=compute_corr(A,s_z)
end
figure
plot(dist, correlator, '-.')

end
function [dist,average]=compute_average(A,loc_op)
average=[];
dist=[];
N = length(A);
for d =1:N
    
    dist(end+1)=d;
    average(end+1) =one_point_corr(A,loc_op,d);
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
function corr_xx=one_point_corr(A,loc_op, point_1)
N=length(A);
A_x =A;

A_x{point_1} =scon({A_x{point_1},loc_op},{[-1,-2,1],[1,-3]});
%A_x{point_2} =scon({A_x{point_2},loc_op},{[-1,-2,1],[1,-3]});

left_x = compute_mix_left(A_x,A,{1});
rigth_x{N} =eye(1);
rigth_x = compute_mix_rigth(A_x,A,rigth_x);
corr_xx = scon({left_x{N+1},rigth_x{N}},{[1,2],[1,2]});
end


