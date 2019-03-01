% Generates a full grid for the evaluation of the exchange-correlation
% integrals in DFT. the function returns the Cartesian coordinates of all
% grid points and their respective weights.
%
% Input:
%    atoms            element numbers, e.g. [6 8] for carbon monoxide
%    xyz_a0           Kx3 array of cartesian coordinates, in bohr radii
%    nRadialPoints    number of radial grid points (typically 50-75)
%    nAngularsPoints  number of angular grid points (possible values:
%                        50, 110, 194, 302, 434, 590, 770)
% Output:
%    grid: structure containing all information about the molecular grid
%      .xyz          Cartesian coordinates for grid points, one per row,
%                      in units of Bohr radii
%      .weights      weights (product of radial and angular weights and r^2)

function grid = molecular_grid(atoms,xyz_a0,nRadialPoints,nAngularPoints)

grid.xyz = [];
grid.weights = [];

for iAtom = 1:numel(atoms)

  % Calculate atomic grid point coordinates and weights
  grid_ = atomic_grid(iAtom,xyz_a0,atoms,nRadialPoints,nAngularPoints);
  
  % Calculate partition weights and total weights for all grid points
  partweights = partitionweights(grid_.xyz,xyz_a0,atoms,iAtom);
  
  % Add atomic grid to molecular grid
  grid.xyz = [grid.xyz; grid_.xyz];
  grid.weights = [grid.weights; partweights.*grid_.weights];
  
end
