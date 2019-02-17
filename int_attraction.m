% Calculates the MxM matrix of electron-nuclear attraction energy integrals.
% 
% Input:
%   atoms       list of element numbers (array with K elements)
%                   e.g. [6 8] for CO 
%   xyz_a0      K×3 array of Cartesian coordinates of nuclei, in bohr
%   basis       basis information, as obtained by buildbasis
% 
% Output:
%   Vne         M×M matrix of attraction energy integrals, in hartrees

function Vne = int_attraction(atoms,xyz_a0,basis)

M = numel(basis);       % number of basis functions in basis
Vne = zeros(M);       % initializing Vne as an MxM matrix of zeros


for C =1:numel(atoms)       % iterates over each atom
    Rc = xyz_a0(C,:);           % nuclear position vector of current atom
    Zc = atoms(C);              % nuclear charge of current atom
    
    for mu = 1:M  % iterates over basis functions
        basis1 = basis(mu);
        A = basis1.A;
        a = basis1.a;
        numPrimA = numel(basis1.alpha);
        
        for nu = 1:M    % iterates over all basis functions again to allow
                        % matrix elements to be calculated
            basis2 = basis(nu);
            B = basis2.A;       % coordinates of basis func 2
            b = basis2.a;       % cartesian exponents of basis func 2
            numPrimB = numel(basis2.alpha);     % number of primitives
            integral = 0;

            for k = 1:numPrimA  % iterating over primitives in basis func 1
                                        
                alpha = basis1.alpha(k); 
                d_k = basis1.d(k);          
                N_k = basis1.N(k);
                 
                for l=1:numPrimB   % iterating over prims in basis func 2
                    
                    beta = basis2.alpha(l);
                    N_l = basis2.N(l);
                    d_l = basis2.d(l);
                    
                    % Defining some variables
                    K_AB = exp(((-alpha*beta)/(alpha+beta))*norm(A-B)^2);
                    p = alpha+beta;
                    P = (alpha.*A + beta.*B)/(alpha + beta);
                    T = p*norm(P-Rc)^2;
                    
                    % Calling the vertical recursion relation function to
                    % evaluate each primitive, then linearly combining them
                    % to evaluate the integral for each matrx element.
                    primitive = VRR(a,b,0,p,P,K_AB,A,B,Rc,T); 
                    integral = integral + d_k*d_l*N_k*N_l*primitive;
                end
            end
            Vne(mu,nu) = Vne(mu,nu) - Zc*integral; % eveluates matrix element
        end
    end
end
end
  
                


