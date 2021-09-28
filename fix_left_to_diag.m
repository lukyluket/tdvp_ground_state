function [A, left, rigth] = fix_left_to_diag(A, left, rigth)
N = length(A);

for n = 1: N
    
    if n ==1
[Ul, el] = eig( 0.5*(left{n} + left{n}'));
         norm = trace(el);
         left{1} = el /norm ;
         A{n} = scon({A{n}, Ul},{[1,-2,-3],[1, -1]});
%left{1} =eye(size(left{1}))/size(left{1},1) ;
new_left =  scon( {  A{n},    conj(A{n}),   left{n}},...
             {[1, -1, 2],   [3 , -2 , 2],     [1, 3]});
%         
        err_gf = abs(new_left - left{n+1}/norm);
     max_err = max( err_gf(1:end));
        left{n+1} = new_left;
    else
        [Ul, el] = eig( 0.5*(left{n} + left{n}'));
        A{n-1} = scon({A{n-1}, Ul'},{[-1,1,-3],[-2,1]});
        
        A{n} = scon({A{n}, Ul},{[1,-2,-3],[1, -1]});
        
        left{n} = scon( {  A{n-1},    conj(A{n-1}),   left{n-1}},...
            {[1, -1, 2],   [3 , -2 , 2],     [1, 3]});
        
        new_left =  scon( {  A{n},    conj(A{n}),   left{n}},...
            {[1, -1, 2],   [3 , -2 , 2],     [1, 3]});
        
%        err_gf = abs(new_left - left{n+1}/norm);
 %       max_err = max( err_gf(1:end));
        left{n+1} = new_left;
    end
    
end
[Ul, el] = eig( 0.5*(left{N+1} + left{N+1}'));
A{N} = scon({A{N}, Ul'},{[-1,1,-3],[-2,1]});
left{N+1} = scon( {  A{N},    conj(A{N}),   left{N}},...
    {[1, -1, 2],   [3 , -2 , 2],     [1, 3]});
rigth{N} = Ul.'*rigth{N}*conj(Ul);
end