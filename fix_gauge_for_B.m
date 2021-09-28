function Vtr = fix_gauge_for_B(A, rigth)
 N =length(A);
    for n = N : (-1) : 1
        R{n}= scon ({conj(A{n}) , sqrt(rigth{n})},...
            {[-1,1,-2],   [-3,1] });
        s1 = size(A{n},1);
        sd = size(A{n},3);
        s2 = size(rigth{n},1);
        R{n} = reshape( R{n}, s1 ,sd*s2);
        [Ur, er, Vr] =svd (R{n});
        Vrd=Vr';
        st =sd*s2-s1;
        if st > 0
        Vtr{n} = reshape(Vrd( end-st+1:end, :),st ,sd ,s2);
        Vtr{n} = permute(Vtr{n}, [1,3, 2]);
        else 
            Vtr{n} = zeros(s1,s2,sd);
        % aVtr = Vr(:,end-st+1:end);
        %aVtr = reshape( aVtr , sd, s2, st);
        %aVtr = permute(aVtr, [3,2,1]) ;
        %qui ci andrebbe un complesso coniugato se lo volessi usare come

        %B pero' poiche lo uso come B^dagger non lo metto cosi' lo tolgo
        %ovunque
        
    end
end
