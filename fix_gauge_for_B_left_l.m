function Vtr = fix_gauge_for_B_left_l(A, left,n)
 N =length(A);
    %for n = N : (-1) : 1
        R{n}= scon ({A{n} , diag(sqrt(diag(left{n})))},...
            {[1,-1,-2],   [-3,1] });
        s1 = size(A{n},2);
        sd = size(A{n},3);
        s2 = size(A{n},1);
        R{n} = reshape( R{n}, s1 ,sd*s2);
        [Ur, er, Vr] =svd (R{n});
        %[Ur, er, Vr] =svdecon(R{n});
        Vrd=Vr';
       % st =sd*s2-s1;
        st =size(Vr,2)-s1;
       if st > 0
        Vtr{n} = reshape(Vr( :,end-st+1:end),sd ,s2,st);
        Vtr{n} = permute(Vtr{n}, [2,3,1]);
        else 
            Vtr{n} = zeros(s2,s1,sd);
        % aVtr = Vr(:,end-st+1:end);
        %aVtr = reshape( aVtr , sd, s2, st);
        %aVtr = permute(aVtr, [3,2,1]) ;
        %qui ci andrebbe un complesso coniugato se lo volessi usare come

        %B pero' poiche lo uso come B^dagger non lo metto cosi' lo tolgo
        %ovunque
        
    %end
        end
end
