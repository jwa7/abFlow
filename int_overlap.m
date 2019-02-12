% Calculates the M×M matrix of overlap integrals. 
%
% Input:
%       basis       basis information, as obtained by buildbasis
% Output:
%       S           M×M matrix of overlap integrals

function S = int_overlap(basis)
    
    M = numel(basis);           % total number of basis funcs in basis
    S = zeros(M);               % initializing S as all-zero MxM matrix
    fac2 = @(n) prod(n:-2:1);   % evaluates the double factorial of n

    for mu=1:M     % iterates over basis functions in basis
        
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
            S_mu_nu = 0;
                                        
            for k=1:numPrimA        % iterating over mu's primitives
            
            alpha_k = basis(mu).alpha(k);   % mu's k-th radial exponent
            d_k = basis(mu).d(k);   % mu's k-th contraction coefficient
            N_k = basis(mu).N(k);  % mu's k-th normalization constant
                                                
                for l=1:numPrimB        % iterating over nu's primitives
                    
                    beta_l = basis(nu).alpha(l);    % nu's l-th radial exponent
                    d_l = basis(nu).d(l);   % nu's l-th contraction coefficient
                    N_l = basis(nu).N(l);  % nu's l-th normalization constant

                    p_kl = alpha_k + beta_l;
                    K_AB = exp(-alpha_k*beta_l*norm(A-B)^2/(alpha_k+beta_l));
                    Integral_3D = K_AB*(pi/p_kl)^(3/2); % initializing the 3D integral 
                                                        % with the prefactors
                    
                    for omega=1:3   % iterating over components x, y, z
                        P_omega = (alpha_k.*A(omega) + beta_l.*B(omega))./(alpha_k + beta_l);
                        I_omega = 0;    % initializing the 1D integral 
                        for i=0:(a(omega)+b(omega))/2
                            q = 2*i;
                            f_q = 0;
                            for j=max(0,q-a(omega)):min(q,b(omega))
                                f_q = f_q + nchoosek(a(omega),q-j)*nchoosek(b(omega),j)*(P_omega-A(omega))^(a(omega)-q+j)*(P_omega-B(omega))^(b(omega)-j);
                            end
                            I_omega = I_omega + f_q*fac2(q-1)/(2*p_kl)^i;
                        end 
                        Integral_3D = Integral_3D * I_omega;    % the 3D intergal is a cumulative 
                                                                % product of the three 1D components
                    end
                    S_mu_nu = S_mu_nu + d_k*d_l * N_k*N_l * Integral_3D;
                end
                
            end
            S(mu, nu) = S_mu_nu;
        end
        
    end
    
end
