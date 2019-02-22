% Calculates the total nuclear repulsion energy in the molecule.
% 
% Input:
%   atoms       list of element numbers (array with K elements) 
%                     e.g. [6 8] for CO 
%   xyz_a0      K×3 array of Cartesian coordinates of nuclei, in bohr
%
% Output:
%   Vnn         total nuclear repulsion energy, in hartrees

function Vnn = nucnucrepulsion(atoms,xyz_a0)

N = numel(atoms);
Vnn = 0;            % initializing the cumulative sum of the nucnucrepulsion

for i = 1:(N-1)     % iterating over nuclei
    for j = i+1:N       % iterating over nuclei again ensuring no overcounting
        Vnn = Vnn + (atoms(i)*atoms(j) / norm(xyz_a0(i,:) - xyz_a0(j,:)));
    end 
end
end

