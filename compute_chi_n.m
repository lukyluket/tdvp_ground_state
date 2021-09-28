function chi = compute_chi_n(N, d, max_chi)
if length(d)==1
chi(1) = 1;
for n = 2 : N/2
chi(n) = min (d*chi(n-1), max_chi);
end

chi(N+1) =1;
for n= N : (-1) :N/2 +1
    
  chi(n) = min (chi(n+1)*d, max_chi);  
end

else
 for n = 2 : N/2+1
     chi(1) = 1;
chi(n) = min (d(n)*chi(n-1), max_chi);
end

chi(N+1) =1;
for n= N : (-1) :N/2 +1
    
  chi(n) = min (chi(n+1)*d(n), max_chi);  
end
end

end