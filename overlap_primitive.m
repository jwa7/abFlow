% Calculates the overlap integral between two unnormalized Gaussian
% primitives.
%
% Input:
%       a, b            Cartesian exponent vectors [ax ay az] and [bx by bz]
%       alpha, beta     radial exponents for the two primitives, inverse bohr squared
%       A,B             center vectors [Ax Ay Az] and [Bx By Bz] in bohr
%
% Output
%       Sp              primitive overlap integral

function Sp = overlap_primitive(a,b,alpha,beta,A,B)

fac2 = @(n) prod(n:-2:1);   % evaluates the double factorial of n

p_kl = alpha + beta;
K_AB = exp(-alpha*beta*norm(A-B)^2/(alpha+beta));
Sp = K_AB*(pi/p_kl)^(3/2); % initializing the 3D integral
                                    % with the prefactors.

for omega=1:3   % iterating over components x, y, z
    
    P_omega = (alpha.*A(omega) + beta.*B(omega))./(alpha + beta);
    I_omega = 0;    % initializing the 1D integral
    
    for i=0:(a(omega)+b(omega))/2
        
        q = 2*i;
        f_q = 0;
        
        for j=max(0,q-a(omega)):min(q,b(omega))
            f_q = f_q + nchoosek(a(omega),q-j)*nchoosek(b(omega),j)*(P_omega-A(omega))^(a(omega)-q+j)*(P_omega-B(omega))^(b(omega)-j);
        end
        
        I_omega = I_omega + f_q*fac2(q-1)/(2*p_kl)^i;
    end
    
    Sp = Sp * I_omega;      % the 3D integral is a cumulative product of the 
                            % three 1D components and the prefactor.
end
end
