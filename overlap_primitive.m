%a) Calculates the overlap integral between two unnormalized Gaussian
%primitives.

%Input:
%a, b = Cartesian exponent vectors [ax ay az] and [bx by bz]
%alpha, beta = radial exponents for the two primitives, inverse bohr squared
%A,B = center vectors [Ax Ay Az] and [Bx By Bz] in bohr

%Output = Sp = primitive overlap integral

function Sp = overlap_primitive(a,b,alpha,beta,A,B)

df = @(n) prod(n:-2:1); %calculates double factorial of n

%Terms for calculating 3D overlap integral that do not 
p = alpha + beta; 
P = (alpha.*A + beta.*B)/(p);
K_AB = exp(((-alpha*beta)/(alpha+beta))*norm(A - B)^2);

Iw = [0 0 0];  %Placeholder

for w = 1:3  %given a = [ax ay az], w loops through ax ay then az, and calculates the 3D overlap integral for each.
    for i=0:(a(w)+b(w)/2) %This loop sums over i=0 to (a(w) + b(w))/2 for calculating Iw, the 1D integral
        k = 2*i;
        fk = 0;
        for j = max(0,k-a(w)) : min(k,b(w))    %This loop calculates fk_temp and sums them into fk, which is used to calculate the 1D integral Iw.    
            binom_coeff1 = factorial(a(w)) / (factorial(k-j)*factorial(a(w) - (k-j))); %Calculates 1st binomial coefficient term
            binom_coeff2 = factorial(b(w)) / (factorial(j)*factorial(b(w)-j)); %Calculates 2nd binomial coefficient term
            term3 = (P(w) - A(w))^(a(w) - k + j); 
            term4 = (P(w) - B(w))^(b(w) - j); 
            fk_temp = binom_coeff1 * binom_coeff2 * term3 * term4; %a single fk term
            fk = fk + fk_temp; %Sum of fk terms 
        end
        Iw_temp = fk*df(2*i-1)/((2*p).^i); %A single 1D integral
        Iw(w) = Iw(w) + Iw_temp; %Sum of the 1D integrals stored according to the appropriate coordinate (ax, ay, or az)
    end
end
Sp = (pi/p)^(3/2) * K_AB * Iw(1) * Iw(2) * Iw(3); %3D overlap integral.



end


