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
%   motest HF           restricted Hartree-Fock SCF (RHF)
%   motest bfeval       basis function evaluation
%   motest grid         integration grid for DFT
%   motest xc           exchange-correlation energy and integrals
%   motest DFT          restricted Kohn-Sham SCF (RKS)
%
%   There are 17 test cases: 1-11 are RHF, 12-17 are RKS

function motest(testname,testid)

% Test error thresholds
%-----------------------------------------------------------------------------
thresholds.Sp = 1e-6; % primitive overlap integrals
thresholds.S = 1e-7; % overlap integrals
thresholds.T = 1e-7; % kinetic-energy integrals
thresholds.BoysF = 1e-10; % Boys function
thresholds.Vne = 1e-7; % attraction integrals
thresholds.ERI = 1e-7; % electron-electron repulsion integrals
thresholds.Vee = 1e-6; % electron-electron repulsion matrix elements
thresholds.Vnn = 1e-7; % nuclear repulsion energy
thresholds.epsilon = 1e-7; % orbital energies
thresholds.FCSCe = 1e-7; % F*C - S*C*epsilon 
thresholds.Cnorms = 1e-7; % norms of MO coefficent vectors
thresholds.Ctheta = 1e-7; % angle between occ. MO subspaces
thresholds.P = 1e-6; % density matrix elements
thresholds.E0 = 1e-6; % electronic energy
thresholds.Etot = 1e-6; % total clamped-nuclei energy
thresholds.Exc = 1e-6; % exchange-correlation energy
thresholds.Vxc = 1e-7; % exchange-correlation integrals
thresholds.bfeval = 1e-8; % basis function evaluation
thresholds.rhoInt = 1e-3; % density integral
%-----------------------------------------------------------------------------

clc

alltests = {'basisdef','basisread','buildbasis',...
  'Sp','S','T','boys','Vne','ERI','Vnn',...
  'HF','bfeval','grid','xc','DFT'};

% Check for test data file
if ~exist('motest.mat','file')
  error('  File motest.mat not found!');
end

% Load test data file
load('motest','testdata','overlapprimitivetest');
for tt = 1:numel(testdata)
  testout(tt) = testdata(tt).results;
end

if nargin<1 || isempty(testname)
  testname = alltests;
end 
if ischar(testname)
  testname = {testname};
end

if nargin<2
  testid = 1:numel(testdata);
end
if nargin==2 && ischar(testid)
  testid = str2double(testid);
  if isnan(testid)
    error('Please specify a single number as the second argument.');
  end
end

if any(testid<1) || any(testid>numel(testdata))
  error('The second argument must be a single number between 1 and %d.',numel(testdata));
end

for iTest = 1:numel(testname)
  
  fprintf('==== motest %s ===========================================\n',testname{iTest});
  switch lower(testname{iTest})
    case 'basisdef'
      pass = motest_basisfiles();
    case 'basisread'
      pass = motest_basisread();
    case 'buildbasis'
      pass = motest_buildbasis(testid);
    case 'sp'
      pass = motest_overlap_primitive();
    case 's'
      pass = motest_overlap(testid);
    case 't'
      pass = motest_kinenergy(testid);
    case 'boys'
      pass = motest_boysF();
    case 'vne'
      pass = motest_attraction(testid);
    case 'eri'
      pass = motest_repulsion(testid);
    case 'vnn'
      pass = motest_nucrep(testid);
    case 'hf'
      pass = motest_hf(testid);
    case 'bfeval'
      pass = motest_bfeval();
    case 'grid'
      pass = motest_molgrid(testid);
    case 'xc'
      pass = motest_xc(testid);
    case 'dft'
      pass = motest_dft(testid);
    otherwise
      error('Test ''%s'' is not known.',testname{iTest});
  end
  
  if pass
    fprintf('  motest %s passed!\n',testname{iTest});
  else
    fprintf('  motest %s FAILED!\n',testname{iTest});
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

b62exp = [7.8683 1.8813 0.5442];
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
      break
    end
    if abs(basis(p).A-basisref(p).A)>1e-3
      ok = false;
      fprintf('    Basis function %d is centered on the wrong atom position.\n',p);
      break
    end
    if any(basis(p).a~=basis(p).a)
      ok = false;
      fprintf('    Basis function %d has incorrect cartesian exponents.\n',p);
      break
    end
    if any(abs(basis(p).alpha-basisref(p).alpha)./basisref(p).alpha>1e-3)
      ok = false;
      fprintf('    Basis function %d has incorrect radial exponents.\n',p);
      break
    end
    if any(abs(basis(p).d-basisref(p).d)>1e-3)
      ok = false;
      fprintf('    Basis function %d has incorrect contraction coefficients.\n',p);
      break
    end
    if any(abs(basis(p).N-basisref(p).N)./basisref(p).N>1e-5)
      ok = false;
      fprintf('    Basis function %d has incorrect primitive normalization constants.\n',p);
      break
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

fprintf('   Sp threshold:   %e\n',thresholds.Sp);

% Run over all pair-wise combinations of s, p, and d functions
for i1 = 1:10
  for i2 = 1:10
    a = spd(i1,:);
    b = spd(i2,:);
    Sref = overlapprimitivetest.S(i1,i2);
    S = overlap_primitive(a,b,alpha,beta,A,B);
    relerr = abs((S - Sref)/Sref);
    fprintf(' <%3s|%3s>   error %e',spd_str{i1},spd_str{i2},relerr);
    if relerr > thresholds.Sp
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

fprintf('  S error threshold:     %e\n',thresholds.S);

for t = testid
  fprintf('  Test ID %2d: ',t);
  bdef = basisread(testdata(t).basisset);
  basis = buildbasis(testdata(t).atoms,testdata(t).xyz_a0,bdef);
  
  S = int_overlap(basis);
  Sref = testout(t).S;
  S_err = maxabsdiff(S,Sref);
  
  fprintf(' max error %e',S_err);
  if S_err > thresholds.S
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end
end

end

%===============================================================================
function ok = motest_kinenergy(testid)

ok = true;

if ~exist('int_kinenergy.m','file')
  ok = false;
  disp('   File int_kinenergy.m is missing.');
  return
end

fprintf('  T error threshold: %e Eh\n',thresholds.T);

for t = testid
  fprintf('  Test ID %2d: ',t);
  bdef = basisread(testdata(t).basisset);
  basis = buildbasis(testdata(t).atoms,testdata(t).xyz_a0,bdef);
  
  T = int_kinenergy(basis);
  Tref = testout(t).T;
  T_err = maxabsdiff(T,Tref);
  
  fprintf(' error %e Eh',T_err);
  if T_err > thresholds.T
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

fprintf('  error threshold:     %e\n',thresholds.BoysF);

ok = true;
for i = 1:numel(m)
  val_ = boysF(m(i),T(i));
  abserr = abs(val_-val(i));
  fprintf('  m=%0.2f,T=%0.2f: error %e',m(i),T(i),abserr);
  if abserr > thresholds.BoysF
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

% Check for presence of correct file
if ~exist('int_attraction.m','file')
  ok = false;
  disp('   File int_attraction.m is missing.');
  return
end

% Test attraction matrix elements for all given testcases
fprintf('  Vne error threshold: %e Eh\n',thresholds.Vne);
for i = testid
  fprintf('  Test ID %2d:   ',i);
  bdef = basisread(testdata(i).basisset);
  basis = buildbasis(testdata(i).atoms,testdata(i).xyz_a0,bdef);
  Vne = int_attraction(testdata(i).atoms,testdata(i).xyz_a0,basis);
  
  Vne_ref = testout(i).Vne;
  Vne_err = maxabsdiff(Vne,Vne_ref);
  
  fprintf(' error %e Eh',Vne_err);
  if Vne_err > thresholds.Vne
    ok = false;
    fprintf(' -- FAIL\n');
  else
    fprintf(' -- pass\n');
  end
end

end

%===============================================================================
function ok = motest_repulsion(testid)

ok = true;

% Check for presence of correct files
files = {'int_repulsion.m','eri_primitive_fast.m'};
for iFile = 1:numel(files)
  if ~exist(files{iFile},'file')
    ok = false;
    fprintf('   File %s is missing.\n',files{iFile});
    return
  end
end

fprintf('  ERI error threshold: %e Eh\n',thresholds.ERI);

for i = testid
  fprintf('  Test ID %2d:   ',i);
  bdef = basisread(testdata(i).basisset);
  basis = buildbasis(testdata(i).atoms,testdata(i).xyz_a0,bdef);
  
  ERI = int_repulsion(basis);
  ERIref = testout(i).ERI;
  ERI_err = maxabsdiff(ERI,ERIref);
  
  fprintf(' error %e Eh',ERI_err);
  if ERI_err > thresholds.ERI
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

fprintf('  Vnn error threshold: %e Eh\n',thresholds.Vnn);

for i = testid
  fprintf('  Test ID %2d:   ',i);
  atoms = testdata(i).atoms;
  xyz_a0 = testdata(i).xyz_a0;
  
  Vnn_ = nucnucrepulsion(atoms,xyz_a0);
  Vnn_ref = nucrepulsionenergy(atoms,xyz_a0);  
  Vnn_err = abs(Vnn_ - Vnn_ref);
  
  fprintf(' error %e Eh',Vnn_err);
  if Vnn_err > thresholds.Vnn
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

ok = true;

xyz_a0 = [0.1 0.5 0.3];
bf.A = [-0.2 -0.5 0.4];
bf.alpha = [1 2 3]/10;
bf.d = [0.4 0.3 0.5];
bf.N = [1.1 1.2 1.3];

a = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 1 0 1; 0 2 0; 0 1 1; 0 0 2];

v_ref = [...
   1.150374213527717, ...
   0.345112264058315, ...
   1.150374213527717, ...
  -0.115037421352772, ...
   0.103533679217495, ...
   0.345112264058315, ...
  -0.034511226405832, ...
   1.150374213527717, ...
  -0.115037421352772, ...
   0.011503742135277];

fprintf('  error threshold:   %e\n',thresholds.bfeval);

str = {'s','px','py','pz','dxx','dxy','dxz','dyy','dyz','dzz'};
for k = 1:size(a,1)
  bf.a = a(k,:);
  v(k) = eval_bf(bf,xyz_a0);
  bf_err = abs(v(k)-v_ref(k));
  if bf_err > thresholds.bfeval
    ok = false;
    fprintf('   %3s   error %e -- FAIL\n',str{k},bf_err);
  else
    fprintf('   %3s   error %e -- pass\n',str{k},bf_err);
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
  nGrid_ = numel(grid.weights);
  nGrid = nRadialPoints*nAngularPoints*numel(atoms);
  if nGrid_~=nGrid
    ok = false;
    fprintf(' -- FAIL\n');
  else
    fprintf(' -- pass\n');
  end
  
end

end

%===============================================================================
function ok = motest_xc(testid)

ok = true;

if ~exist('int_xc.m','file')
  ok = false;
  disp('   File int_repulsion.m is missing.');
  return
end

fprintf('  Vxc error threshold:   %e\n',thresholds.Vxc);

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
  Vxc_err = maxabsdiff(Vxc,Vxc_ref);
  
  fprintf('     Vxc    error %e',Vxc_err);
  if Vxc_err > thresholds.Vxc
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end
  
  % Test Exc, exchange-correlation energy
  Exc_ref = testout(t).Exc;
  Exc_err = abs(Exc-Exc_ref);
  
  fprintf('     Exc    error %e',Exc_err);
  if Exc_err > thresholds.Exc
    fprintf(' -- FAIL\n');
    ok = false;
  else
    fprintf(' -- pass\n');
  end

  % Test rhoInt, integral of electron density
  rhoInt_ref = sum(atoms)-testdata(t).charge;
  rhoInt_err = abs(rhoInt-rhoInt_ref);
  
  fprintf('     rhoInt error %e',rhoInt_err);
  if rhoInt_err > thresholds.rhoInt
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
  fields = {'basis','S','T','Vne','J','K','ERI','epsilon','C','P','E0','Etot'};
  if isDFT
    fields = [fields {'Vxc','Exc','rhoInt'}];
  end
  for f = 1:numel(fields)
    if ~isfield(out,fields{f})
      ok = false;
      fprintf('     out.%s is missing -- FAIL\n',fields{f});
    end
  end
  
  % Compare electron-electron repulsion matrix, Vee = J - K
  if isfield(out,'J') && isfield(out,'K')
    Vee = out.J - out.K;
    Vee_err = maxabsdiff(Vee,testout(t).Vee);
    if Vee_err > thresholds.Vee
      ok = false;
      fprintf('     J-K       error %e Eh -- FAIL\n',Vee_err);
    else
      fprintf('     J-K       error %e Eh -- pass\n',Vee_err);
    end
  end
  
  % Compare exchange-correlation matrix, Vxc
  if isDFT
    if isfield(out,'Vxc')
      Vxc_err = maxabsdiff(out.Vxc,testout(t).Vxc);
      if Vxc_err > thresholds.Vxc
        ok = false;
        fprintf('     Vxc       error %e Eh -- FAIL\n',Vxc_err);
      else
        fprintf('     Vxc       error %e Eh -- pass\n',Vxc_err);
      end
    end
  end
  
  % Check size of orbital energy vector, epsilon
  if isfield(out,'epsilon')
    eps_correctSize = numel(out.epsilon)==numel(testout(t).epsilon);
    if ~eps_correctSize
      ok = false;
      fprintf('     epsilon   error: incorrect size -- FAIL\n');
    else
      fprintf('     epsilon   correct size          -- pass\n');
    end
  end
  
  % Check sorting of orbital energies, epsilon
  if isfield(out,'epsilon') && eps_correctSize
    if any(diff(out.epsilon(:))<0)
      ok = false;
      fprintf('     epsilon   sorted incorrectly    -- FAIL\n');
    else
      fprintf('     epsilon   sorted correctly      -- pass\n');
    end
  end
  
  % Compare orbital energies, epsilon
  if isfield(out,'epsilon') && eps_correctSize
    eps_err = maxabsdiff(out.epsilon,testout(t).epsilon);
    if eps_err > thresholds.epsilon
      ok = false;
      fprintf('     epsilon   error %e Eh -- FAIL\n',eps_err);
    else
      fprintf('     epsilon   error %e Eh -- pass\n',eps_err);
    end
  end
  
  % Test whether C has the correct size
  if isfield(out,'C')
    C_correctSize = true;
    if ~isequal(size(out.C),size(testout(t).C))
      C_correctSize = false;
      ok = false;
      fprintf('     C         error: wrong size   -- FAIL\n');
    end
  end
  
  % Test whether FC == SCe is satisfied
  if isfield(out,'C') && C_correctSize
    F = testout(t).T + testout(t).Vne + testout(t).Vee;
    if isDFT
      F = F + testout(t).Vxc;
    end
    S = testout(t).S;
    C = out.C;
    epsilon = diag(testout(t).epsilon);
    maxabserr = maxabsdiff(F*C,S*C*epsilon);
    if maxabserr > thresholds.FCSCe
      ok = false;
      fprintf('     FC-SCe    error %e Eh -- FAIL\n',maxabserr);
    else
      fprintf('     FC-SCe    error %e Eh -- pass\n',maxabserr);
    end
  end
  
  % Check C normalization
  if isfield(out,'C') && C_correctSize
    C = out.C;
    S = out.S;
    Cnorms = diag(C'*S*C);
    maxerr_Cnorms = max(abs(Cnorms-1));
    if maxerr_Cnorms > thresholds.Cnorms
      ok = false;
      fprintf('     C norms   error %e    -- FAIL\n',maxerr_Cnorms);
    else
      fprintf('     C norms   error %e    -- pass\n',maxerr_Cnorms);
    end
  end
  
  % Compare subspace spanned by the occupied MO coefficient vectors
  if isfield(out,'C') && C_correctSize
    nElectrons = sum(testdata(t).atoms)-testdata(t).charge;
    occ = 1:nElectrons/2;
    C_occ = out.C(:,occ);
    C_occ_ref = testout(t).C(:,occ);
    theta = subspace(C_occ,C_occ_ref); % angle between occupied subspaces, ideally zero
    if theta > thresholds.Ctheta
      ok = false;
      fprintf('     C space   error %e    -- FAIL\n',theta);
    else
      fprintf('     C space   error %e    -- pass\n',theta);
    end
  end
    
  % Compare density matrix, P
  if isfield(out,'P')
    P_err = maxabsdiff(out.P,testout(t).P);
    if P_err > thresholds.P
      ok = false;
      fprintf('     P         error %e    -- FAIL\n',P_err);
    else
      fprintf('     P         error %e    -- pass\n',P_err);
    end
  end
  
  % Compare electronic energy, E0
  if isfield(out,'E0')
    E0_err = abs(out.E0-testout(t).E0);
    if E0_err > thresholds.E0
      ok = false;
      fprintf('     E0        error %e Eh -- FAIL\n',E0_err);
    else
      fprintf('     E0        error %e Eh -- pass\n',E0_err);
    end
  end
  
  % Compare total energy, Etot
  if isfield(out,'Etot')
    Etot_err = abs(out.Etot-testout(t).Etot);
    if Etot_err > thresholds.Etot
      ok = false;
      fprintf('     Etot      error %e Eh -- FAIL\n',Etot_err);
    else
      fprintf('     Etot      error %e Eh -- pass\n',Etot_err);
    end
  end
  
  % Check whether density integral is correct
  if isDFT &&isfield(out,'rhoInt')
    nElectrons = sum(testdata(t).atoms) - testdata(t).charge;
    rhoInt_err = abs(out.rhoInt-nElectrons);
    if rhoInt_err > thresholds.rhoInt
      ok = false;
      fprintf('     rhoInt    error %e    -- FAIL\n',rhoInt_err);
    else
      fprintf('     rhoInt    error %e    -- pass\n',rhoInt_err);
    end
  end
  
  % Compare exchange-correlation energy, Exc
  if isDFT && isfield(out,'Exc')
    Exc_err = maxabsdiff(out.Exc,testout(t).Exc);
    if Exc_err > thresholds.Exc
      ok = false;
      fprintf('     Exc       error %e Eh -- FAIL\n',Exc_err);
    else
      fprintf('     Exc       error %e Eh -- pass\n',Exc_err);
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

% For DFT, add information on integration grid
isDFT = strcmp(t.method,'RKS');
if isDFT
  str = [str sprintf('|%dang|%drad',t.nAngularPoints,t.nRadialPoints)];
end

end

%===============================================================================
function mad = maxabsdiff(A,B)
mad = max(abs(A(:)-B(:)));
end

%===============================================================================
function  Vnn = nucrepulsionenergy(Z,xyz_a0)
Vnn = 0;
nAtoms = numel(Z);
for a = 1:nAtoms
  for b = a+1:nAtoms
    A = xyz_a0(a,:);
    B = xyz_a0(b,:);
    Vnn = Vnn + Z(a)*Z(b)/norm(A-B);
  end
end
end
