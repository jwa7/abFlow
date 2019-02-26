function val = eval_bf(basisfun,xyz_a0)
a = basisfun.a;
alpha = basisfun.alpha;
N = basisfun.N;
A = basisfun.A;
coeffs = basisfun.d;
evalbf = zeros(size(xyz_a0,1),1);

for i=1:size(xyz_a0,1)
   
    for j=1:length(alpha);
        evalbf(i) = evalbf(i) + N(j)*coeffs(j).*exp(-alpha(j)*(sum((xyz_a0(i,:)-A).^2)));
    end
    evalbf(i) = evalbf(i)*((xyz_a0(i,1)-A(1)).^a(1) .* (xyz_a0(i,2)-A(2)).^a(2) .* (xyz_a0(i,3)-A(3)).^a(3));   

end
val = evalbf;
end
