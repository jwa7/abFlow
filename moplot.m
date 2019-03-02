% Plots a specific molecular orbital.
% 
% Input:
%     atoms       list of element numbers (1×K array); e.g. [6 8] for CO
%     xyz_a0      coordinates of the nuclei in the molecule, in bohr
%     out         output structure as returned from the mocalc function
%     iMO         index of MO to plot (iMO = 1 for lowest-energy MO, 
%                 etc. in ascending-energy order)
%     level       contour level for the isosurface, in bohr-based units a0?(3/2)
%     
% Output: 
%     A 3D plot of the MO, including positions of the atoms.

function moplot(atoms,xyz_a0,out,iMO,level)

nAtoms = numel(atoms);                  % number of atoms in molecule.
dimensions = max(max(xyz_a0)) + 5;      % defining the plotting dimensions.
C = out.C;                              % the MO coefficient matrix.
basis = out.basis;                      % the basis functions for the molecule.
M = numel(basis);                       % number of basis functions.
points = 100;                           % number of grid points.

x = linspace(-dimensions, dimensions, points);
y = linspace(-dimensions, dimensions, points);
z = linspace(-dimensions, dimensions, points);
[X,Y,Z] = meshgrid(x,y,z);

bf = zeros(points, points, points);

for mu = 1:M
    for x1 = 1:points
        for y1 = 1:points
            for z1 = 1:points
                point = [x(x1) y(y1) z(z1)];
                bf(x1,y1,z1) = bf(x1,y1,z1) + C(mu,iMO)*eval_bf(basis(mu), point);
            end
        end
    end
end

figure
hold on

for i = 1:nAtoms
    scatter3(xyz_a0(i,1), xyz_a0(i,2), xyz_a0(i,3), 5^3*atoms(i), 'filled');
end

isosurface(X,Y,Z,bf,level);
isosurface(X,Y,Z,bf,-level);

axis equal;
title(sprintf('Molecular Orbital Diagram for MO number %d at isolevel %d', iMO, level));
xlabel('x / a_0');
ylabel('y / a_0');
zlabel('z / a_0');
lighting gouraud
camlight('headlight')
alpha(0.4);
view(3)

end