function I3D = integral_3D(alpha_k,beta_l,A,B,a,b)

fac2 = @(n) prod(n:-2:1);   % evaluates the double factorial of n

p_kl = alpha_k + beta_l;
K_AB = exp(-alpha_k*beta_l*norm(A-B)^2/(alpha_k+beta_l));
Integral_3D = K_AB*(pi/p_kl)^(3/2); % initializing the 3D integral
                                    % with the prefactors.

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
    
    Integral_3D = Integral_3D * I_omega;    % the 3D integral is a cumulative
                                            % product of the three 1D components
end
I3D = Integral_3D;
end