function weights = partitionweights(xyzGrid,xyzAtoms,Atoms,atomIdx)

% xyzGrid: cartesian coordinates of all grid points (nGridx3 array)
% xyzAtom: cartesian coordaintes of all atoms (nAtomsx3 array)
% atomIdx: index of center for which the partition weights should be
%  calculated

% Calculates the partition weights accoring to the 'fuzzy Voronoi'
% approach of A. D. Becke. Equation numbers in the following refer to
%   Becke, J. Chem. Phys. 88, 2547 (1988), http://dx.doi.org/10.1063/1.454033

nAtoms = size(xyzAtoms,1);
nGrid = size(xyzGrid,1);

p = @(mu) 1.5*mu - 0.5*mu.^3; % base polynomial, Eq.(19)
s = @(mu) (1 - p(p(p(mu))))/2; % 3rd-order cutoff profile, Eqs.(20) and (21)

% Calculate distances of each grid point to all atomic centers
R = zeros(nGrid,nAtoms);
for iPt = 1:nGrid
  for iAtom = 1:nAtoms
    R(iPt,iAtom) = norm(xyzGrid(iPt,:)-xyzAtoms(iAtom,:));
  end
end

% Calculate cell functions P for all grid points
for iAtom = 1:nAtoms
  P_ = ones(nGrid,1);
  for jAtom = 1:nAtoms
    if jAtom == iAtom, continue; end
    RIJ = norm(xyzAtoms(iAtom,:)-xyzAtoms(jAtom,:));
    muIJ = (R(:,iAtom) - R(:,jAtom)) / RIJ; % Eq.(11)
    chi = atomicradius(Atoms(iAtom))/atomicradius(Atoms(jAtom)); % Eq.(A4)
    aIJ = (1-chi^2)/(4*chi); % Eqs.(A5) and (A6) combined
    aIJ = max(min(aIJ,0.5),-0.5); % Eq.(A3)
    P_ = P_.*s(muIJ + aIJ*(1-muIJ.^2)); % Eqs. (13), (A1), (A2)
  end
  P(:,iAtom) = P_;
end

% Calculate partition weights
weights = P(:,atomIdx)./sum(P,2);  % Eq.(22)

function R = atomicradius(iElements)

% List of covalent radii (Cordero, 2008; based on 37k crystal structures)
% DOI: http://dx.doi.org/10.1039/B801115J
R = [0.31 0.28 1.28 0.96 0.84 (0.76+0.73+0.69)/3 0.71 0.66 0.57 0.58]; % H-He, Li-Ne
% (Becke uses Bragg-Slater radii (Slater, 1964))

R = R(iElements);

bohr = 0.52917721067; % Bohr radius, in Angstrom
R = R / bohr; % Angstrom -> bohr
