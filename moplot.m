% Plots a specific molecular orbital.
% 
% Input:
%     atoms       list of element numbers (1×K array); e.g. [6 8] for CO
%     xyz_a0      coordinates of the nuclei in the molecule, in bohr
%     out
%     iMO
%     level
% Output: 
%     A 3D plot of the MO, including positions of the atoms.

function moplot(atoms,xyz_a0,out,iMO,level)

nAtoms = numel(atoms);                  % number of atoms in molecule.
dimensions = max(max(xyz_a0)) + 5;      % defining the plotting dimensions.
M = numel(out.basis);                   % number of basis functions.
nucSize = 10;                           % size of atom on plot.
limits = [-dimensions dimensions];      % dimensional limits for plot.

for nuc = 1:nAtoms
    hold on;
    plot3(xyz_a0(nuc,1), xyz_a0(nuc,2), xyz_a0(nuc,2))
%     , "o", "MarkerSize", nucSize, "FaceAlpha", level);
    set(gca, 'Xlim', limits, 'Ylim', limits, 'Zlim', limits);
    view(3);
end

x = -dimensions : dimensions;
y = -dimensions : dimensions;
z = -dimensions : dimensions;
[X,Y,Z] = meshgrid(x,y,z);
bf_coords = [X(:), Y(:), Z(:)];

bf = zeros(1,numel(bf_coords));

for basis = 1:M
   bf = bf + MO_eval_bf(out.basis(basis),bf_coords, out.C(basis,iMO)); 
end


end