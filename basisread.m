% Reads in a basis set definition file (in Gaussian94 format)
% and returns a structure containing all basis set details.
%
%  Input:
%    BasisSetName   basis set name, e.g. '6-31G'
%
%  Output:
%    BasisSet         cell array with one entry per element
%       BasisSet{iElement} is a structure array containing the fields
%         .shelltype     type of shell ('S', 'P', 'SP', 'D', etc)
%         .coeffs        contraction coefficients
%         .exponents     radial exponents of primitives
%
%  The basis set files are expected to be in a subfolder called
%  basissets. The file names are expected to have the extension
%  .basis. For instance, the 6-31G basis should be in the file
%  basissets/6-31G.basis

function BasisSet = basisread(BasisSetName)

basissetfilename = ['basissets/' BasisSetName '.basis'];

% Read all lines from file into a single cell array (one cell per line)
f = fopen(basissetfilename);
if f<0
  error('Could not open file ''%s''.',basissetfilename);
end
allLines = textscan(f,'%s','Delimiter','\n','Whitespace','');
allLines = allLines{1};
fclose(f);

% Remove comment lines (starting with !) and empty lines
for iLine = 1:numel(allLines)
  rmv(iLine) = isempty(allLines{iLine}) || (allLines{iLine}(1)=='!');
end
allLines(rmv) = [];

% Identify lines delimiting atoms (****)
for iLine = 1:numel(allLines)
  delimLines(iLine) = (allLines{iLine}(1)=='*');
end
delimLines = find(delimLines);
nEntries = numel(delimLines)-1; % number of entries (i.e. atoms) in the file

for iEntry = 1:nEntries
  
  % Extract lines for this atom
  lines = allLines(delimLines(iEntry)+1:delimLines(iEntry+1)-1);
  
  % Get element number
  ElementNo = elementsymb2no(lines{1}(1:2));
  
  % Determine number of shells
  shellDelim = 0;
  for iLine = 2:numel(lines)
    shellDelim(iLine) = lines{iLine}(1)~=' ';
  end
  shellDelim(end+1) = true;
  shellDelim = find(shellDelim);
  nShells = numel(shellDelim)-1;

  % Read in each shell in turn
  for iShell = 1:nShells
    % Get shell type (S, P, SP, D, etc) from shell header line
    shellHeader = lines{shellDelim(iShell)};
    newShell.shelltype = strtrim(shellHeader(1:2));
    % Get exponents and contraction coefficients
    stri = char(lines(shellDelim(iShell)+1:shellDelim(iShell+1)-1));
    numbers = str2num(stri);  %#ok
    exponents = numbers(:,1).';
    coeffs = numbers(:,2:end).';
    % Store exponents and contraction coefficients in output structure
    newShell.exponents = exponents;
    newShell.coeffs = coeffs;
    BasisSet{ElementNo}(iShell) = newShell;
  end
end


function elno = elementsymb2no(symbol)

if numel(symbol)==1, symbol(2) = ' '; end

ElementSymbols = [...
  'H He'...
  'LiBeB C N O F Ne'...
  'NaMgAlSiP S ClAr'...
  'K CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr'...
  'RbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI Xe'...
  'CsBaLaCePrNdPmSmEuGdTbDyHoerTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn'...
  'FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRgCnNhFlMcLvTsOg'];

idx = strfind(ElementSymbols,symbol);
elno = (idx+1)/2;
