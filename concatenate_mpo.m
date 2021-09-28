
function [c_mpo_left, c_mpo_rigth, O] = concatenate_mpo(l_c_mpo_left,l_c_mpo_rigth, l_O );
 
KK = length(l_O);
offset =0;
for k =1 :KK
last_el = size(l_O{k}{1},1);
c_mpo_left(1,offset+1:offset+last_el) =l_c_mpo_left{k};
c_mpo_rigth(offset+1:offset+last_el,1) =l_c_mpo_rigth{k};
for n =1 :length(l_O{k})
O{n}(offset +1 :offset+last_el,offset +1:offset+last_el, :,: ) = l_O{k}{n};
end
offset = offset + last_el;
end

