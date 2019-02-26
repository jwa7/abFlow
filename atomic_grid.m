% grid = atomic_grid(iAtom,xyz_a0,atoms,nRadialPoints,nAngularPoints)
%
% Generates a full 3D spherical grid (product of 1D radial and
% 2D angular grid) centered at Center and radially scaled for
% the element given in Element.
%
% It returns the Cartesian coordinates of the grid points and
% their respective weights, taking into account the partition weighing
% within the given molecule.
% 
% Input:
%    iAtom            index of atom for which the grid should be calculated
%    atom             element numbers, e.g. [6 8] for carbon monoxide
%    xyz_a0           Kx3 array of cartesian coordinates, in bohr radii
%    nRadialPoints    number of radial grid points (typically 50-75)
%    nAngularsPoints  number of angular grid points (possible values:
%                        50, 110, 194, 302, 434, 590, 770)
%
% Output:
%    grid: structure containing all information about the atomic grid
%      .xyz          Cartesian coordinates for grid points, one per row,
%                      in units of Bohr radii
%      .weights      weights (product of radial and angular weights and r^2)

function grid = atomic_grid(iAtom,xyz_a0,atoms,nRadialPoints,nAngularPoints)

Element = atoms(iAtom);
Center = xyz_a0(iAtom,:);

% Define atomic radii
%--------------------------------------------------------
% List of covalent radii (Cordero, 2008; based on 37k crystal structures)
% DOI: http://dx.doi.org/10.1039/B801115J
R(1:2) = [0.31 0.28]; % H-He, Angstrom
R(3:10) =  [1.28 0.96 0.84 (0.76+0.73+0.69)/3 0.71 0.66 0.57 0.58]; % Li-Ne, Angstrom
R = R/0.52917721067;  % Angstrom -> bohr
rm = R(Element)/2; % scaling factor for grid

% Generate radial grid: Chebyshev-Gauss (first kind)
%--------------------------------------------------------
xi = cos((2*(1:nRadialPoints).'-1)/(2*nRadialPoints)*pi);
w = ones(nRadialPoints,1)*pi/nRadialPoints;
% mapping: (-1,1) -> (0,inf), including scaling by rm
radialGrid.r = rm*(1+xi)./(1-xi);
radialGrid.weights = rm*2*(sqrt(1-xi.^2)./(1-xi).^2) .* w;

% Generate angular grid: Lebedev
%--------------------------------------------------------
angularGrid = lebedev_grid(nAngularPoints);

% Combine into full product grid
%--------------------------------------------------------
% combine coordinates
x = radialGrid.r * angularGrid.x.';
y = radialGrid.r * angularGrid.y.';
z = radialGrid.r * angularGrid.z.';
% shift to requested center
grid.xyz = [x(:)+Center(1), y(:)+Center(2), z(:)+Center(3)];

% calculate total weights
r2 = x(:).^2 + y(:).^2 + z(:).^2;
grid.weights = radialGrid.weights * angularGrid.w.';
grid.weights = grid.weights(:).*r2;
