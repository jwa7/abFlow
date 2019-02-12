% Calculates the M×M matrix of kinetic-energy integrals.
% 
% Input:
%     basis       basis information, as obtained by buildbasis
% Output:
%     T           M×M matrix of kinetic-energy integrals, in hartrees

function T = int_kinenergy(basis)

M = numel(basis);           % total number of basis funcs in basis
T = zeros(M);               % initializing T as all-zero MxM matrix
fac2 = @(n) prod(n:-2:1);   % evaluates the double factorial of n

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
        T_mu_nu = 0;
        
        for k=1:numPrimA        % iterating over mu's primitives
            
            alpha_k = basis(mu).alpha(k);   % mu's k-th radial exponent
            d_k = basis(mu).d(k);   % mu's k-th contraction coefficient
            N_k = basis(mu).N(k);  % mu's k-th normalization constant
            
            for l=1:numPrimB        % iterating over nu's primitives
                
                beta_l = basis(nu).alpha(l);    % nu's l-th radial exponent
                d_l = basis(nu).d(l);   % nu's l-th contraction coefficient
                N_l = basis(nu).N(l);  % nu's l-th normalization constant
                
                I3D_kinetic = 0;
                
                for omega=1:3   % iterating over components x, y, z
                    I_omega = 0;
                    
                    if b >= 2
                        I_omega = beta_l*(2*b(omega)+1)*integral_3D(alpha_k,beta_l,A,B,a,b) ...
                                  -2*beta_l^2*integral_3D(alpha_k,beta_l,A,B,a,b+2);
                                  -0.5*b(omega)*(b(omega)-1)*integral_3D(alpha_k,beta_l,A,B,a,b-2);
                    end
                    
                    I3D_kinetic = I3D_kinetic + I_omega;
                    
                end
                T_mu_nu = T_mu_nu + d_k*d_l * N_k*N_l * I3D_kinetic;
            end
        end
        T(mu, nu) = T_mu_nu;
    end
    
end
end