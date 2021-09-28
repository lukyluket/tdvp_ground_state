function [c_mpo_left, c_mpo_rigth, O] = build_mpo(Oid,Ox, Oy, Oz,coup_c,N)
sx = [0, 1;...
      1, 0];
sz = [1, 0;...
    0, -1];
sy = [0, 1i;...
     -1i, 0];
 id =eye(2);


    MO(:,:,1)= id;
     MO(:,:,2)= sx;
      MO(:,:,3)= sy;
     MO(:,:,4)= sz;
     
     A(:,:,1) =Oid;
     A(:,:,2) =Ox;
    A(:,:,3) =Oy;
     A(:,:,4) =Oz;
     
     O{1}= scon({A, MO},{[-1,-2,1],[-4,-3,1]});
     for n=2:N
     O{n}= O{1};
     end
     c_mpo_left = zeros(1,size(Oid,1));
     c_mpo_left(1) =coup_c;
     c_mpo_rigth = zeros(size(Oid,2),1);
     c_mpo_rigth(end) =1;

end