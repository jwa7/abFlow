function primitive = VRR(a,b,m,p,P,K_AB,A,B,C,T)

if sum([a b]<0) > 0 % returns 0 if have a negative value
    
    primitive = 0;
    return
    
elseif a==b & a==[0 0 0]    % no recursive reduction needed if a = b = [0 0 0]
    
    primitive = (2*pi*K_AB/p)*boysF(m,T);
    return
    
else
    
    z = [0 0 0];
    
    % Recursively reduces a to [0 0 0]
    red_a = @(w,z)(P(w)-A(w))*VRR(a-z,b,m,p,P,K_AB,A,B,C,T) ...
            + (C(w)-P(w))*VRR(a-z,b,m+1,p,P,K_AB,A,B,C,T) ... 
            +((a(w)-1)/(2*p))*(VRR(a-z*2,b,m,p,P,K_AB,A,B,C,T) ... 
            - VRR(a-z*2,b,m+1,p,P,K_AB,A,B,C,T)) ...
            +(b(w)/(2*p))*(VRR(a-z,b-z,m,p,P,K_AB,A,B,C,T) ... 
            - VRR(a-z,b-z,m+1,p,P,K_AB,A,B,C,T));
    
    % Recursively reduces b to [0 0 0]
    red_b = @(w,z)(P(w)-B(w))*VRR(a,b-z,m,p,P,K_AB,A,B,C,T) ... 
            + (C(w)-P(w))*VRR(a,b-z,m+1,p,P,K_AB,A,B,C,T) ...
            + (b(w)-1)/2*p * (VRR(a,b-z*2,m,p,P,K_AB,A,B,C,T) ... 
            - VRR(a,b-z*2,m+1,p,P,K_AB,A,B,C,T)) ...
            + a(w)/2*p * (VRR(a-z,b-z,m,p,P,K_AB,A,B,C,T) ... 
            - VRR(a-z,b-z,m+1,p,P,K_AB,A,B,C,T));
    
    for w=1:3 % If a(w) doesnt =0, inputs a into red_a.
        
        if a(w)>0
            
            z(w) = 1;
            primitive = red_a(w,z);
            return
            
        end
    end
    
    for w=1:3 % If b(w) doesnt equal zero, inputs a into red_b.
        
        if b(w)>0
            
            z(w) = 1;
            primitive = red_b(w,z);
            return
            
        end
    end
end