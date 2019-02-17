function recursion = VRR(a,b,m,p,P,K_AB,A,B,C,T)

if sum([a b]<0) > 0 %returns 0 if have a negative value
    recursion = 0;
    return
end

if a==b & a==[0 0 0] %calculates if a = b = [000]
    recursion = (2*pi*K_AB/p)*boysF(m,T);
    return
end  

z = [0 0 0];

red_a = @(w,z)(P(w)-A(w))*VRR(a-z,b,m,p,P,K_AB,A,B,C,T) + (C(w)-P(w))*VRR(a-z,b,m+1,p,P,K_AB,A,B,C,T) ... %Reduces a to [0 0 0]
    +((a(w)-1)/(2*p))*(VRR(a-z*2,b,m,p,P,K_AB,A,B,C,T) - VRR(a-z*2,b,m+1,p,P,K_AB,A,B,C,T)) ...
    +(b(w)/(2*p))*(VRR(a-z,b-z,m,p,P,K_AB,A,B,C,T) - VRR(a-z,b-z,m+1,p,P,K_AB,A,B,C,T));

red_b = @(w,z)(P(w)-A(w))*VRR(a,b-z,m,p,P,K_AB,A,B,C,T) + (C(w)-P(w))*VRR(a,b-z,m+1,p,P,K_AB,A,B,C,T) ... %Reduces b to [0 0 0]
    +(a(w)-1)/2*p * (VRR(a,b-z*2,m,p,P,K_AB,A,B,C,T) - VRR(a,b-z*2,m+1,p,P,K_AB,A,B,C,T)) ...
    +b(w)/2*p * (VRR(a-z,b-z,m,p,P,K_AB,A,B,C,T) - VRR(a-z,b-z,m+1,p,P,K_AB,A,B,C,T));

for w=1:3 %If a(w) doesnt =0, inputs a into red_a.
    if a(w)>0
        z(w) = 1;
        recursion = red_a(w,z);
        return
    end
end

for w=1:3 %If b(w) doesnt =0, inputs a into red_b.
    if b(w)>0
        z(w) = 1;
        recursion = red_b(w,z);
        return
    end
end
end