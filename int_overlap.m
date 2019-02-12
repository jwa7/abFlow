% Calculates the M×M matrix of overlap integrals. 
%
% Input:
%       basis       basis information, as obtained by buildbasis
% Output:
%       S           M×M matrix of overlap integrals

function S = int_overlap(basis)

M = numel(basis);           % total number of basis funcs in basis
S = zeros(M);               % initializing S as all-zero MxM matrix

for mu=1:M     % iterates over all basis functions in basis
    
    A = basis(mu).A;
    a = basis(mu).a;
    numPrimA = numel(basis(mu).alpha);  % number of primitives in
                                        % the mu-th basis function.
    
    for nu=1:M      % iterates over basis funcs in basis again to
                    % compare each func with every other one
        
        B = basis(nu).A;
        b = basis(nu).a;
        numPrimB = numel(basis(nu).alpha);  % number of primitives in
                                            % nu-th basis function.
        
        for k=1:numPrimA        % iterating over mu's primitives
            
            alpha = basis(mu).alpha(k);   % mu's k-th radial exponent
            d_k = basis(mu).d(k);   % mu's k-th contraction coefficient
            N_k = basis(mu).N(k);  % mu's k-th normalization constant
            
            for l=1:numPrimB        % iterating over nu's primitives
                
                beta = basis(nu).alpha(l);    % nu's l-th radial exponent
                d_l = basis(nu).d(l);   % nu's l-th contraction coefficient
                N_l = basis(nu).N(l);  % nu's l-th normalization constant
                
                Sp = overlap_primitive(a,b,alpha,beta,A,B); 
                
                S(mu,nu) = S(mu,nu) + d_k*d_l * N_k*N_l * Sp; % matrix element 
                                                              % is a cumulative sum
            end
        end
    end
end
end
