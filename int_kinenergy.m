% Calculates the M×M matrix of kinetic-energy integrals.
% 
% Input:
%     basis       basis information, as obtained by buildbasis
% Output:
%     T           M×M matrix of kinetic-energy integrals, in hartrees

function T = int_kinenergy(basis)

M = numel(basis);           % total number of basis funcs in basis
T = zeros(M,M);               % initializing T as all-zero MxM matrix

for mu=1:M     % iterates over all basis functions in basis
    basis1 = basis(mu);
    A = basis1.A;
    a = basis1.a;
    numPrimA = numel(basis1.alpha);  % number of primitives in
                                        % the mu-th basis function.
    
    for nu=1:M      % iterates over basis funcs in basis again to
                    % compare each func with every other one
        basis2 = basis(nu);
        B = basis2.A;
        b = basis2.a;
        numPrimB = numel(basis2.alpha);  % number of primitives in
                                            % nu-th basis function.
        
        for k=1:numPrimA        % iterating over mu's primitives
            
            alpha = basis1.alpha(k);   % mu's k-th radial exponent
            d_k = basis1.d(k);   % mu's k-th contraction coefficient
            N_k = basis1.N(k);  % mu's k-th normalization constant
            
            for l=1:numPrimB        % iterating over nu's primitives
                
                beta = basis2.alpha(l);    % nu's l-th radial exponent
                d_l = basis2.d(l);   % nu's l-th contraction coefficient
                N_l = basis2.N(l);  % nu's l-th normalization constant
                
                Iw = [0 0 0];
                for w=1:3   % iterating over components x, y, z
                   
                    prim2plus = b;
                    prim2minus = b;
                    prim2plus(w) = b(w)+2; %Adds 2 to the w-th value in b
                    prim2minus(w) = b(w)-2; %Subtracts 2 from the w-th value in b
                    
                    integral1 = overlap_primitive(a,b,alpha,beta,A,B);  %Calculates overlap integral [k,a|l,b]                                    
                    integral2 = overlap_primitive(a,prim2plus,alpha,beta,A,B); %Calculates overlap integral [k,a|l,b+2]
                    
                    if b < 0 %Sets [k,a|,b-2] = 0 if b<0
                        integral3 = 0; 
                    else
                        integral3 = overlap_primitive(a,prim2minus,alpha,beta,A,B); %Calculates overlap integral [k,a|l,b-2]
                    end  
                    
                    Iw(w) = beta*(2*b(w)+1)*integral1 - 2*beta^2 * integral2 - 1/2*b(w)*(b(w)-1)*integral3;                  
                end
                T(mu,nu) = T(mu,nu) + d_k*d_l * N_k*N_l * (Iw(1)+Iw(2)+Iw(3));
            end
        end
    end
end
end