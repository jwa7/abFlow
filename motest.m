% motest  Test function for all aspects of RHF and RKS implementations
%
%   motest basisdef     basis set definition files
%   motest basisread    reading of basis set definition files
%   motest buildbasis   building basis for molecules
%   motest Sp           primitive overlap integrals
%   motest S            overlap integrals
%   motest T            kinetic-energy integrals
%   motest boys         Boys function
%   motest Vne          electron-nuclear attraction energy integrals
%   motest ERI          electron-electron repulsion energy integrals
%   motest Vnn          nuclear-nuclear repulsion energy
%   motest HF           Hartree-Fock SCF (RHF)
%   motest bfeval       basis function evaluation
%   motest grid         integration grid for DFT
%   motest xc           exchange-correlation energy and integrals
%   motest DFT          Kohn-Sham SCF (RKS)
%
%   There are 17 test cases: 1-11 are RHF, 12-17 are RKS

function motest(testname,testid)

clc

alltests = {'basisdef','basisread','buildbasis',...
  'Sp','S','T','boys','Vne','ERI','Vnn',...
  'HF','bfeval','grid','xc','DFT'};

% Check for test data files
if ~exist('motest.mat','file')
  error('  File motest.mat not found!');
end

% Load test data
load('motest','testdata','overlapprimitivetest');
for tt = 1:numel(testdata)
  testout(tt) = testdata(tt).results;
end

if nargin<1 || isempty(testname), testname = alltests; end 
if nargin<2, testid = 1:numel(testdata); end
if ischar(testname), testname = {testname}; end

for iTest = 1:numel(testname)
  
  fprintf('==== motest %s ===========================================\n',testname{iTest});
  switch testname{iTest}
    case 'basisdef'
      pass = motest_basisfiles();
    case 'basisread'
      pass = motest_basisread();
    case 'buildbasis'
      pass = motest_buildbasis(testid);
    case 'Sp'
      pass = motest_overlap_primitive();
    case 'S'
      pass = motest_overlap(testid);
    case 'T'
      pass = motest_kinenergy(testid);
    case 'boys'
      pass = motest_boysF();
    case 'Vne'
      pass = motest_attraction(testid);
    case 'ERI'
      pass = motest_repulsion(testid);
    case 'Vnn'
      pass = motest_nucrep(testid);
    case {'HF','hf'}
      pass = motest_hf(testid);
    case 'bfeval'
      pass = motest_bfeval();
    case 'grid'
      pass = motest_molgrid(testid);
    case {'XC','xc'}
      pass = motest_xc(testid);
    case {'DFT','dft'}
      pass = motest_dft(testid);
    otherwise
      error('Test ''%s'' is not known.',testname{iTest});
  end
  
  if pass
    fprintf('motest %s passed!\n',testname{iTest});
  else
    fprintf('motest %s FAILED!\n',testname{iTest});
  end
    
end

%===============================================================================
function ok = motest_basisfiles
  
basissets = {'STO-3G','6-31G','6-311G','cc-pVDZ'};

ok = true;

if ~exist('basissets','dir')
  ok = false;
  fprintf('   Folder basissets/ is missing.\n');  
end

for b = 1:numel(basissets)
  filename = ['basissets/' basissets{b} '.basis'];
  if ~exist(filename,'file')
    ok = false;
    fprintf('   File %s is missing.\n',filename);
  end
end

end

%===============================================================================
function ok = motest_basisread

ok = true;
if ~exist('basisread.m','file')
  ok = false;
  disp('   File basisread.m is missing.');
  return
end

b = basisread('6-31G');

b62exp = [7.8683    1.8813    0.5442];
okbasis = numel(b)==10 && numel(b{6})==3 && all(abs(b{6}(2).exponents-b62exp)<1e-4);

if ~okbasis
  disp('    Basis set definition file basissets/6-31G.basis is not correct.');
end

ok = ok && okbasis;

end

%===============================================================================
function ok = motest_buildbasis(testid)

ok = true;
if ~exist('buildbasis.m','file')
  ok = false;
  disp('   File buildbasis.m is missing.');
  return
end

for t = testid
  teststr = testsummary(testdata(t));
  fprintf('  Test %2d: %s \n',t,teststr);
  
  bdef = basisread(testdata(t).basisset);
  atoms = testdata(t).atoms;
  xyz_a0 = testdata(t).xyz_a0;
  basis = buildbasis(atoms,xyz_a0,bdef);
  basisref = testout(t).basis;
  
  M = numel(basis);
  Mref = numel(basisref);
  
  % Verify that we have the correct number of basis functions
  if M~=Mref
    ok = false;
    fprintf('   Number of basis functions for %s is not correct: %d given, %d expected.\n',...
      testdata(t).molecule,M,Mref);
    break
  end
  
  % Verify that all the fields are present
  fieldnames = {'atom','A','a','alpha','d','N'};
  for ifn = 1:numel(fieldnames)
    if ~isfield(basis,fieldnames{ifn})
      ok = false;
      error('The basis function list missed the field %s.',fieldnames{ifn});
    end
  end
  
  % Verify that all fields have the correct type and size
  for p = 1:M
    for ifn = 1:numel(fieldnames)
      fn = fieldnames{ifn};
      if ~isa(basis(p).(fn),class(basisref(p).(fn)))
        error('The field basis(%d).%s has the wrong type.',p,fn);
      end
      if ~isequal(size(basis(p).(fn)),size(basisref(p).(fn)))
        error('The field basis(%d).%s has the wrong size.',p,fn);
      end
    end
  end
  
  % Verify that we have the correct ordering of basis functions
  alphasums = @(ba) arrayfun(@(x)sum(x.alpha),ba);
  sortkeys = @(ba) [cat(1,ba.atom) cat(1,ba.A) cat(1,ba.a) alphasums(ba).'];  
  skref = sortkeys(basisref);
  sk = sortkeys(basis);
  
  % Sort test and reference basis functions
  [~,idx] = sortrows(sk);
  [~,idxref] = sortrows(skref);
  if any(idx~=idxref)
    ok = false;
    fprintf('    Basis functions are not in the required order (by atoms, shells, cartesian exponents).\n');
    break
  end
    
  % Verify that basis function are correct
  for p = 1:M
    if basis(p).atom~=basisref(p).atom
      ok = false;
      fprintf('    Basis function %d is centered on the wrong atom type.\n',p);
      break;
    end
    if abs(basis(p).A-basisref(p).A)>1e-3
      ok = false;
      fprintf('    Basis function %d is centered on the wrong atom position.\n',p);
      break;
    end
    if any(basis(p).a~=basis(p).a)
      ok = false;
      fprintf('    Basis function %d has incorrect cartesian exponents.\n',p);
      break;
    end
    if any(abs(basis(p).alpha-basisref(p).alpha)./basisref(p).alpha>1e-3)
      ok = false;
      fprintf('    Basis function %d has incorrect radial exponents.\n',p);
      break;
    end
    if any(abs(basis(p).d-basisref(p).d)>1e-3)
      ok = false;
      fprintf('    Basis function %d has incorrect contraction coefficients.\n',p);
      break;
    end
    if any(abs(basis(p).N-basisref(p).N)./basisref(p).N>1e-5)
      ok = false;
      fprintf('    Basis function %d has incorrect primitive normalization constants.\n',p);
      break;
    end
  end
  
  if ~ok, break; end
  
end

end

%===============================================================================
function ok = motest_overlap_primitive

ok = true;
if ~exist('overlap_primitive.m','file')
  ok = false;
  disp('   File overlap_primitive.m is missing.');
  return
end
  
A = overlapprimitivetest.A;
B = overlapprimitivetest.B;
alpha = overlapprimitivetest.alpha;
beta = overlapprimitivetest.beta;
spd = overlapprimitivetest.ab;

spd_str = {'s','px','py','pz','dx2','dxy','dxz','dy2','dyz','dzz'};

threshold = 1e-6;
fprintf('   threshold:              %e\n',threshold);

% Run over all pair-wise combinations of s, p, and d functions
for i1 = 1:10
  for i2 = 1:10
    a = spd(i1,:);
    b = spd(i2,:);
    Sref = overlapprimitivetest.S(i1,i2);
    S = overlap_primitive(a,b,alpha,beta,A,B);
    relerr = abs((S - Sref)/Sref);
    fprintf(' <%3s|%3s>   error %e',spd_str{i1},spd_str{i2},relerr);
    if relerr > threshold
      ok = false;
      fprintf(' -- FAIL\n');
    else
      fprintf(' -- pass\n');
    end
  end
end

end

%===============================================================================
function ok = motest_overlap(testid)

ok = true;
if ~exist('int_overlap.m','file')
  ok = false;
  disp('   File int_overlap.m is missing.');
  return
end

threshold = 1e-7;
fprintf('  error threshold:             %e\n',threshold);

for t = testid
  fprintf('  Test ID %2d: ',t);
  bdef = basisread(testdata(t).basisset);
  basis = buildbasis(testdata(t).atoms,testdata(t).xyz_a0,bdef);
  S = int_overlap(basis);
  
  Sref = testout(t).S;
  abserr = (S(:)-Sref(:));
  maxabserr(t) = max(abserr(:));
  ok = all(abserr(:)<threshold);
  
  fprintf(' max error %e',maxabserr(t));
  if maxabserr(t)>threshold
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end
end
fprintf('  max abs error: %e\n',max(maxabserr));

end

%===============================================================================
function ok = motest_kinenergy(testid)

ok = true;
if ~exist('int_kinenergy.m','file')
  ok = false;
  disp('   File int_kinenergy.m is missing.');
  return
end

threshold = 1e-7;
fprintf('  error threshold:   %e Eh\n',threshold);

for t = testid
  fprintf('  Test ID %2d: ',t);
  bdef = basisread(testdata(t).basisset);
  basis = buildbasis(testdata(t).atoms,testdata(t).xyz_a0,bdef);
  T = int_kinenergy(basis);
  
  Tref = testout(t).T;
  abserr = (T(:)-Tref(:));
  maxabserr = max(abserr(:));
  ok = all(abserr(:)<threshold);
  
  fprintf(' error %e Eh',maxabserr);
  if maxabserr>threshold
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end
end

end

%===============================================================================
function ok = motest_boysF

if ~exist('boysF.m','file')
  ok = false;
  disp('   File boysF.m is missing.');
  return
end

mTval(1,:) = [1 2 0.115702180856173];
mTval(2,:) = [2 0.5 0.140750536825913];
mTval(3,:) = [0 3.2 0.489762207784512];
mTval(4,:) = [0 0 1];
mTval(5,:) = [1 0 1/3];
mTval(6,:) = [2 0 0.2];

m = mTval(:,1);
T = mTval(:,2);
val = mTval(:,3);

threshold = 1e-10;
ok = true;
for i = 1:numel(m)
  val_ = boysF(m(i),T(i));
  abserr = abs(val_-val(i));
  fprintf('  m=%g,T=%g: error %e',m(i),T(i),abserr);
  if abserr>threshold
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end
end

end

%===============================================================================
function ok = motest_attraction(testid)

ok = true;
if ~exist('int_attraction.m','file')
  ok = false;
  disp('   File int_attraction.m is missing.');
  return
end

threshold = 1e-7;
fprintf(' error threshold:    %e\n',threshold);

for i = testid
  fprintf('  Test ID %2d: ',i);
  bdef = basisread(testdata(i).basisset);
  basis = buildbasis(testdata(i).atoms,testdata(i).xyz_a0,bdef);
  Vne = int_attraction(testdata(i).atoms,testdata(i).xyz_a0,basis);
  
  Vne_ref = testout(i).Vne;
  abserr = abs(Vne(:)-Vne_ref(:));
  maxabserr = max(abserr(:));
  ok = all(abserr(:)<threshold);
  
  fprintf(' error %e',maxabserr);
  if maxabserr>threshold
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end
end

end

%===============================================================================
function ok = motest_repulsion(testid)

if ~exist('int_repulsion.m','file')
  ok = false;
  disp('   File int_repulsion.m is missing.');
  return
end

ok = true;
if ~exist('eri_primitive_fast.m','file')
  ok = false;
  disp('   File eri_primitive_fast.m is missing.');
  return
end

threshold = 1e-7;
fprintf('  error threshold:   %e\n',threshold);

for i = testid
  fprintf('  Test ID %2d: ',i);
  bdef = basisread(testdata(i).basisset);
  basis = buildbasis(testdata(i).atoms,testdata(i).xyz_a0,bdef);
  ERI = int_repulsion(basis);
  ERI_ref = testout(i).ERI;
  abserr = (ERI(:)-ERI_ref(:));
  maxabserr = max(abserr(:));
  ok = all(abserr(:)<threshold);
  
  fprintf(' error %e',maxabserr);
  if maxabserr>threshold
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end
end

end

%===============================================================================
function ok = motest_nucrep(testid)

ok = true;
if ~exist('nucnucrepulsion.m','file')
  ok = false;
  disp('   File nucnucrepulsion.m is missing.');
  return
end

threshold = 1e-7;
fprintf('  error threshold:   %e Eh\n',threshold);

for i = testid
  fprintf('  Test ID %2d: ',i);
  atoms = testdata(i).atoms;
  xyz_a0 = testdata(i).xyz_a0;
  
  Vnn = 0;
  nAtoms = numel(atoms);
  for a = 1:nAtoms
    for b = a+1:nAtoms
      A = xyz_a0(a,:);
      B = xyz_a0(b,:);
      Vnn = Vnn + atoms(a)*atoms(b)/norm(A-B);
    end
  end
  Vnn_ = nucnucrepulsion(atoms,xyz_a0);
  
  abserr = abs(Vnn_ - Vnn);
  
  fprintf(' error %e Eh',abserr);
  if abserr>threshold
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end
end

end

%===============================================================================
function ok = motest_hf(testid)

ok = motest_scf(testid,'HF');
  
end


%===============================================================================
function ok =  motest_bfeval()

xyz_a0 = [0.1 0.5 0.3];
bf.A = [-0.2 -0.5 0.4];
bf.alpha = [1 2 3]/10;
bf.d = [0.4 0.3 0.5];
bf.N = [1.1 1.2 1.3];

a = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 1 0 1; 0 2 0; 0 1 1; 0 0 2];

v_ref = [   1.150374213527717, ...
   0.345112264058315, ...
   1.150374213527717, ...
  -0.115037421352772, ...
   0.103533679217495, ...
   0.345112264058315, ...
  -0.034511226405832, ...
   1.150374213527717, ...
  -0.115037421352772, ...
   0.011503742135277];

threshold = 1e-8;

ok = true;
str = {'s','px','py','pz','dxx','dxy','dxz','dyy','dyz','dzz'};
for k = 1:size(a,1)
  bf.a = a(k,:);
  v(k) = eval_bf(bf,xyz_a0);
  err(k) = abs(v(k)-v_ref(k));
  if err(k)>threshold
    ok = false;
    fprintf('   %3s   error %e -- FAIL\n',str{k},err(k));
  else
    fprintf('   %3s   error %e -- pass\n',str{k},err(k));
  end
end

end

%===============================================================================
function ok = motest_molgrid(testid)
  
filenames = {'molecular_grid.m','atomic_grid.m','partitionweights.m','lebedev_grid.m'};
ok = true;
for iFile  = 1:numel(filenames)
  if ~exist(filenames{iFile},'file')
    ok = false;
    fprintf('   File %s is missing.\n',filenames{iFile});
  end
end

for t = testid
  teststr = testsummary(testdata(t));
  isDFT = strcmp(testdata(t).method,'RKS');
  if ~isDFT
    fprintf('  Test %2d: %s -- skipped\n',t,teststr);
    continue
  else
    fprintf('  Test %2d: %s',t,teststr);
  end
  
  atoms = testdata(t).atoms;
  xyz_a0 = testdata(t).xyz_a0;
  nRadialPoints = testdata(t).nRadialPoints;
  nAngularPoints = testdata(t).nAngularPoints;
  grid = molecular_grid(atoms,xyz_a0,nRadialPoints,nAngularPoints);
  n_ = numel(grid.weights);
  n = nRadialPoints*nAngularPoints*numel(atoms);
  if n_~=n
    ok = false;
    fprintf(' -- FAIL\n');
  else
    fprintf(' -- pass\n');
  end
  
end

end

%===============================================================================
function ok = motest_xc(testid)

if ~exist('int_xc.m','file')
  ok = false;
  disp('   File int_repulsion.m is missing.');
  return
end

ok = true;

threshold = 1e-7;
fprintf('  error threshold:   %e\n',threshold);

for t = testid
  teststr = testsummary(testdata(t));
  isDFT = strcmp(testdata(t).method,'RKS');
  if ~isDFT
    fprintf('  Test %2d: %s -- skipped\n',t,teststr);
    continue
  else
    fprintf('  Test %2d: %s \n',t,teststr);
  end
  
  bdef = basisread(testdata(t).basisset);
  basis = buildbasis(testdata(t).atoms,testdata(t).xyz_a0,bdef);
  atoms = testdata(t).atoms;
  xyz_a0 = testdata(t).xyz_a0;
  Exch = testdata(t).ExchFunctional;
  Corr = testdata(t).CorrFunctional;
  nRadial = testdata(t).nRadialPoints;
  nAngular = testdata(t).nAngularPoints;
  P = testout(t).P;
  
  grid = molecular_grid(atoms,xyz_a0,nRadial,nAngular);
  [Vxc,Exc,rhoInt] = int_xc(basis,P,grid,Exch,Corr);
  
  % Test Vxc, matrix of exchange-correlation integrals
  Vxc_ref = testout(t).Vxc;
  abserr = abs(Vxc(:)-Vxc_ref(:));
  maxabserr = max(abserr(:));
  
  fprintf('     Vxc    error %e',maxabserr);
  if maxabserr>threshold
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end

  % Test Exc, exchange-correlation energy
  Exc_ref = testout(t).Exc;
  abserr = abs(Exc-Exc_ref);
  
  fprintf('     Exc    error %e',abserr);
  if abserr>threshold
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end

  % Test rhoInt, integral of electron density
  threshold = 1e-3;
  rhoInt_ref = sum(atoms)-testdata(t).charge;
  abserr = abs(rhoInt-rhoInt_ref);
  
  fprintf('     rhoInt error %e',abserr);
  if abserr>threshold
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end
  
end

end

%===============================================================================
function ok = motest_dft(testid)

ok = motest_scf(testid,'DFT');
  
end

%===============================================================================
function ok = motest_scf(testid,testmethod)

ok = true;
if ~exist('mocalc.m','file')
  ok = false;
  disp('   File mocalc.m is missing.');
  return
end

runDFT = strcmp(testmethod,'DFT');

for t = testid
  teststr = testsummary(testdata(t));
  isDFT = strcmp(testdata(t).method,'RKS');
  if xor(runDFT,isDFT)
    fprintf('  Test %2d: %s -- skipped\n',t,teststr);
    continue
  else
    fprintf('  Test %2d: %s \n',t,teststr);
  end
  
  % Assemble input and call mocalc()
  atoms = testdata(t).atoms;
  charge = testdata(t).charge;
  xyz_a0 = testdata(t).xyz_a0;
  settings.method = testdata(t).method;
  settings.basisset = testdata(t).basisset;
  settings.tolEnergy = testdata(t).tolEnergy;
  settings.tolDensity = testdata(t).tolDensity;
  settings.nRadialPoints = testdata(t).nRadialPoints;
  settings.nAngularPoints = testdata(t).nAngularPoints;
  settings.ExchFunctional = testdata(t).ExchFunctional;
  settings.CorrFunctional = testdata(t).CorrFunctional;
  out = mocalc(atoms,xyz_a0,charge,settings);
  
  % Check output structure for all required fields
  fields = {'basis','S','T','Vne','Vee','ERI','epsilon','C','P','E0','Etot'};
  if isDFT
    fields = [fields {'Vxc','Exc','rhoInt'}];
  end
  for f = 1:numel(fields)
    if ~isfield(out,fields{f})
      ok = false;
      fprintf('     out.%s is missing\n',fields{f});
    end
  end
  
  % Compare orbital energies, epsilon
  if isfield(out,'epsilon')
    eps_threshold = 1e-7;
    eps_err = max(abs(out.epsilon-testout(t).epsilon));
    if eps_err>eps_threshold
      ok = false;
      fprintf('     epsilon error %e Eh -- FAIL\n',eps_err);
    else
      fprintf('     epsilon error %e Eh -- pass\n',eps_err);
    end
  else
      fprintf('     epsilon missing!              -- FAIL\n');
  end 
    
  % Compare MO coefficient matrix, C
  if ~isfield(out,'C')
    ok = false;
    fprintf('     C       missing!              -- FAIL\n');
  elseif ~isequal(size(out.C),size(testout(t).C))
    ok = false;
    fprintf('     C       error: wrong size   -- FAIL\n');
  else
    nElectrons = sum(testdata(t).atoms)-testdata(t).charge;
    occ = 1:nElectrons/2;
    C_occ = out.C(:,occ);
    C_occ_ref = testout(t).C(:,occ);
    theta = subspace(C_occ,C_occ_ref); % angle between occupied subspaces, ideally zero
    C_threshold = 1e-7;
    if theta > C_threshold
      ok = false;
      fprintf('     C       error %e    -- FAIL\n',theta);
    else
      fprintf('     C       error %e    -- pass\n',theta);
    end
  end
  
  % Compare density matrix, P
  if ~isfield(out,'P')
    fprintf('     P       missing!              -- FAIL\n');
  else
    P_threshold = 1e-6;
    P_err = max(abs(out.P(:)-testout(t).P(:)));
    if P_err > P_threshold
      ok = false;
      fprintf('     P       error %e    -- FAIL\n',P_err);
    else
      fprintf('     P       error %e    -- pass\n',P_err);
    end
  end
  
  % Compare electronic energy, E0
  if isfield(out,'E0')
    E_threshold = 1e-6;
    E0_err = abs(out.E0-testout(t).E0);
    if E0_err>E_threshold
      ok = false;
      fprintf('     E0      error %e Eh -- FAIL\n',E0_err);
    else
      fprintf('     E0      error %e Eh -- pass\n',E0_err);
    end
  else
      fprintf('     E0      missing!              -- FAIL\n');
  end
  
  % Compare total energy, Etot
  if isfield(out,'Etot')
    E_threshold = 1e-6;
    Etot_err = max(abs(out.Etot-testout(t).Etot));
    if Etot_err>E_threshold
      ok = false;
      fprintf('     Etot    error %e Eh -- FAIL\n',Etot_err);
    else
      fprintf('     Etot    error %e Eh -- pass\n',Etot_err);
    end
  else
      fprintf('     Etot     missing!              -- FAIL\n');
  end
  
  % Compare electron-electron repulsion matrix, Vee
  if isfield(out,'Vee')
    Vee_threshold = 1e-6;
    Vee_err = max(abs(out.Vee(:)-testout(t).Vee(:)));
    if Vee_err>Vee_threshold
      ok = false;
      fprintf('     Vee     error %e Eh -- FAIL\n',Vee_err);
    else
      fprintf('     Vee     error %e Eh -- pass\n',Vee_err);
    end
  else
      fprintf('     Vee     missing!              -- FAIL\n');
  end
  
  % Compare exchange-correlation matrix, Vxc
  if isDFT
    if isfield(out,'Vxc')
      Vxc_threshold = 1e-6;
      Vxc_err = max(abs(out.Vxc(:)-testout(t).Vxc(:)));
      if Vxc_err>Vxc_threshold
        ok = false;
        fprintf('     Vxc     error %e Eh -- FAIL\n',Vxc_err);
      else
        fprintf('     Vxc     error %e Eh -- pass\n',Vxc_err);
      end
    else
      fprintf('     Vxc     missing!              -- FAIL\n');
    end
  end
  
  % Check whether density integral is correct
  if isDFT
    if isfield(out,'rhoInt')
      nElectrons = sum(testdata(t).atoms) - testdata(t).charge;
      rhoInt_threshold = 1e-5;
      err = abs(out.rhoInt-nElectrons);
      if err>rhoInt_threshold
        ok = false;
        fprintf('     rhoInt  error %e    -- FAIL\n',err);
      else
        fprintf('     rhoInt  error %e    -- pass\n',err);
      end
    else
      fprintf('     rhoInt   missing!              -- FAIL\n');
    end
  end
  
  % Compare exchange-correlation energy, Exc
  if isDFT
    if isfield(out,'Exc')
      E_threshold = 1e-6;
      Exc_err = max(abs(out.Exc-testout(t).Exc));
      if Exc_err>E_threshold
        ok = false;
        fprintf('     Exc     error %e Eh -- FAIL\n',Exc_err);
      else
        fprintf('     Exc     error %e Eh -- pass\n',Exc_err);
      end
    else
      fprintf('     Vxc     missing!              -- FAIL\n');
    end
  end
  
end % for

end % motest_dft

end % motest
%===============================================================================
%===============================================================================


%===============================================================================
function str = testsummary(t)
str = [t.molecule '|' t.method '|' t.basisset];
if length(str)<16, str(16) = ' '; end
end
