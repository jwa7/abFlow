% ERI = int_repulsion(basis)
%
% Computes the electron-electron repulsion energy integrals for a given
% basis and returns a 4D array with the results.
%
%   Input:
%     basis     1xM list of basis functions with the following fields
%        .atom    element number
%        .A       cartesian coordinates of center (Bohr radii)
%        .a       cartesian exponents [ax ay az]
%        .alpha   list of radial exponents (inverse Bohr radii)
%        .d       list of contraction coefficients
%        .N       list of normalization constants
%
%   Output:
%     ERI       MxMxMxM 4D array of repulsion integrals

function ERI = int_repulsion(basis)

nBasis = numel(basis);
ERI = zeros(nBasis,nBasis,nBasis,nBasis);

% Loop over all unique combinations of 4 basis functions (takes
% symmetry into account)
for mu = 1:nBasis
  for nu = 1:mu
    for kap = 1:mu
      for lam = 1:kap
        
        % Copy needed information from basis function list, for
        % faster access
        A = basis(mu).A;
        B = basis(nu).A;
        C = basis(kap).A;
        D = basis(lam).A;
        alpha = basis(mu).alpha;
        beta =  basis(nu).alpha;
        gamma = basis(kap).alpha;
        delta = basis(lam).alpha;
        
        dNmu = basis(mu).d.*basis(mu).N;
        dNnu = basis(nu).d.*basis(nu).N;
        dNlam = basis(kap).d.*basis(kap).N;
        dNsig = basis(lam).d.*basis(lam).N;
        
        a = basis(mu).a;
        b = basis(nu).a;
        c = basis(kap).a;
        d = basis(lam).a;
        
        % Calculate repulsion integral <mu,nu|lam,sig> by contracting
        % primitive integrals [k,l|n,o]
        eri_ = 0;
        for k = 1:numel(dNmu)
          for l = 1:numel(dNnu)
            for n = 1:numel(dNlam)
              for o = 1:numel(dNsig)
                eri_ = eri_ + ...
                    dNmu(k)*dNnu(l)*dNlam(n)*dNsig(o) * ...
                      eri_primitive_fast(A,B,C,D,...
                           alpha(k),beta(l),gamma(n),delta(o),...
                           a,b,c,d);
              end
            end
          end
        end
        
        % Utilize symmetry properties of two-electron integrals
        ERI(mu,nu,kap,lam) = eri_;
        ERI(nu,mu,kap,lam) = eri_;
        ERI(mu,nu,lam,kap) = eri_;
        ERI(nu,mu,lam,kap) = eri_;
        
        ERI(kap,lam,mu,nu) = eri_;
        ERI(kap,lam,nu,mu) = eri_;
        ERI(lam,kap,mu,nu) = eri_;
        ERI(lam,kap,nu,mu) = eri_;
        
      end
    end
  end
end


