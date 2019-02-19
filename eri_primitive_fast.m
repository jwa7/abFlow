function Integral = eri_primitive_fast(A,B,C,D,alpha,beta,gamma,delta,a,b,c,d)

p = alpha + beta;
q = gamma + delta;

P = (alpha*A + beta*B)/p;
Q = (gamma*C + delta*D)/q;
W = (p*P + q*Q)/(p+q);

AB = A-B;
CD = C-D;

% Pre-calculate [00|00]^(m) auxiliary integrals
T = p*q/(p+q)*norm(P-Q).^2;
m = 0:(sum(a) + sum(b) + sum(c) + sum(d));
KAB = exp(-alpha*beta /p*sum(AB.^2));
KCD = exp(-gamma*delta/q*sum(CD.^2));
ssss_m = 2*pi^(5/2)/p/q/sqrt(p+q)*KAB*KCD*boysF(m,T);

PA = P-A;
WP = W-P;
QC = Q-C;
WQ = W-Q;

w = 1; Cx = eri_coeffs_1D(a(w),b(w),c(w),d(w));
w = 2; Cy = eri_coeffs_1D(a(w),b(w),c(w),d(w));
w = 3; Cz = eri_coeffs_1D(a(w),b(w),c(w),d(w));
Coeffs = conv2(conv2(Cx,Cy),Cz);

Integral = Coeffs * ssss_m.';


  function result = eri_coeffs_1D(a,b,c,d)
    
    if (b>=1)
      % HRR-b: Transfers angular momentum from b to a
      %------------------------------------------------------------------
      % [a,b|c,d]^(m) = [a+1,b-1|c,d]^(m) + AB(i)*[a,b-1|c,d]^(m)
      t1 = eri_coeffs_1D(a+1,b-1,c,d);
      t2 = eri_coeffs_1D(a  ,b-1,c,d);
      result = t1 + AB(w)*[t2 0];
      
    elseif (d>=1)
      % HRR-d: Transfers angular momentum from d to c
      %------------------------------------------------------------------
      % [a,b|c,d]^(m) = [a,b|c+1,d-1]^(m) + CD(i)*[a,b|c,d-1]^(m)
      t1 = eri_coeffs_1D(a,b,c+1,d-1);
      t2 = eri_coeffs_1D(a,b,c  ,d-1);
      result = t1 + CD(w)*[t2 0];
      
    elseif (a>=1)
      % VRR-a: Reduces angular momentum a
      %------------------------------------------------------------------
      % [a,0|c,0]^(m) = PA(i)*[a-1,0|c,0]^(m) + WP(i)*[a-1,0|c,0]^(m+1)
      %                 + (a-1)/(2p)([a-2,0|c,0]^(m) - q/(p+q)*[a-2,0|c,0]^(m+1)])
      %                 + c/2/(p+q)*[a-1,0|c-1,0]^(m+1)
      t = eri_coeffs_1D(a-1,0,c,0);
      result = PA(w)*[t 0] + WP(w)*[0 t];
      if (a>=2)
        t = eri_coeffs_1D(a-2,0,c,0);
        v = (a-1)/(2*p)*([t 0] - q/(p+q)*[0 t]);
        result = result + [v 0];
      end
      if (c>=1)
        t = eri_coeffs_1D(a-1,0,c-1,0);
        v = c/2/(p+q)*[0 t];
        result = result + [v 0];
      end
      
    elseif (c>=1)
      % VRR-c: Reduces angular momentum c
      %------------------------------------------------------------------
      % [a,0|c,0]^(m) = QC(x)*[a,0|c-1,0]^(m) + WQ(x)*[a,0|c-1,0]^(m+1)
      %                 + (c-1)/(2q)([a,0|c-2,0]^(m) - p/(p+q)*[a,0|c-2,0]^(m+1)])
      %                 + a/2/(p+q)*[a-1,0|c-1,0]^(m+1)
      t = eri_coeffs_1D(a,0,c-1,0);
      result = QC(w)*[t 0] + WQ(w)*[0 t];
      if (c>=2)
        t = eri_coeffs_1D(a,0,c-2,0);
        v = (c-1)/(2*q)*([t 0] - p/(p+q)*[0 t]);
        result = result + [v 0];
      end
      if (a>=1)
        t = eri_coeffs_1D(a-1,0,c-1,0);
        v = a/2/(p+q)*[0 t];
        result = result + [v 0];
      end
      
    else
      
      % [0,0|0,0]^(m)
      %------------------------------------------------------------------
      % (m-dependent prefactors are included outside of recursion)
      result = 1;
      
    end
  end



end