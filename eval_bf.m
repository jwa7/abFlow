% Evaluates a basis function at a list of positions in space.
% 
% Input:
%     basisfun    structure with basis function information, 
%                 from list of basis functions.
%     xyz_a0      mx3 list of Cartesian coordinates, one point per row, 
%                 over which to evaluate the basis function.
% Output:
%     val         m-element array, with the values of the basis function.

function val = eval_bf(basisfun,xyz_a0)
a = basisfun.a;
alpha = basisfun.alpha;
N = basisfun.N;
A = basisfun.A;
coeffs = basisfun.d;
evalbf = zeros(size(xyz_a0,1),1);

for i=1:size(xyz_a0,1)
   
    for j=1:length(alpha)
        evalbf(i) = evalbf(i) + N(j)*coeffs(j).*exp(-alpha(j)*(sum((xyz_a0(i,:)-A).^2)));
    end
    evalbf(i) = evalbf(i)*((xyz_a0(i,1)-A(1)).^a(1) .* (xyz_a0(i,2)-A(2)).^a(2) .* (xyz_a0(i,3)-A(3)).^a(3));   

end
val = evalbf;
end
