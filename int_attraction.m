function V_ne = int_attraction(atoms,xyz_a0,basis)

M = numel(basis);
V_ne = zeros(M,M);


for w = 1:numel(atoms);
    C = xyz_a0(w,:);
    
    for mu = 1:M  %basis 1
        basis1 = basis(mu);
        A = basis1.A;
        a = basis1.a;
        numPrimA = numel(basis1.alpha);
        
        for nu = 1:M %basis 2
            basis2 = basis(nu);
            B = basis2.A;
            b = basis2.a;
            numPrimB = numel(basis2.alpha);
            basissum = 0;

            for k = 1:numPrimA %primitive 1     
                alpha = basis1.alpha(k);
                 d_k = basis1.d(k);
                 N_k = basis1.N(k);
                 
                
                for l=1:numPrimB %primitive 2
                beta = basis2.alpha(l);
                N_l = basis2.N(l);
                d_l = basis2.d(l);
                 
                K_AB = exp(((-alpha*beta)/(alpha+beta))*norm(A-B)^2); %constants/variables
                p = alpha+beta;
                P = (alpha.*A + beta.*B)/(alpha + beta);
                T = p*norm(P-C)^2;
                
                
                primitive = VRR(a,b,0,p,P,K_AB,A,B,C,T); %calculates primitive
                basissum = basissum + d_k*d_l*N_k*N_l*primitive; %basis set & sums them
                end
            end
            V_ne(mu,nu) = V_ne(mu,nu) - atoms(w)*basissum; %calculates V_ne & sums
        end
    end
end
end
  
                


