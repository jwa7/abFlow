% int_xc calculates the matrix of exchange-correlation integrals in Vxc, the
% exchange-correlation energy in Exc, and the integral over the electron
% density in rhoInt.

function [Vxc,Exc,rhoInt] = int_xc(basis,P,grid,ExchFunctional,CorrFunctional)
M = length(basis);

for mu=1:M %Evalutes all basis functions at each grid point
    basisfxn(:,mu) = eval_bf(basis(mu),grid.xyz);
end

rhotemp = zeros(length(grid.xyz),1);

for mu = 1:M %Calculates product of two basis functions and density rho
    for nu = 1:M
        rhotemp = rhotemp + basisfxn(:,mu).*basisfxn(:,nu).*P(mu,nu);
    end
end

rho = rhotemp;
rhoInt = sum(grid.weights .* rho);

%Slater exchange
Cx = (3/4)*(3/pi)^(1/3);
epsi_x = -Cx*rho.^(1/3);
Vx = -4/3 * Cx * rho.^(1/3);

%Vosko-Wilk-Nusair correlation
if CorrFunctional == 'VWM3'
    b = 13.0720; %VWN3
    c = 42.7198;
    x0 = -0.4099286;
else 
    b = 3.72744; %VWN5
    c = 12.9352;
    x0 = -0.10498;
end

%Constants/definitions and calculations for VWM Correlation
A = 0.0310907; %Eh
x = (3/(4*pi*rho)).^(1/6);
Q = sqrt(4*c-b^2);
eta = atan(Q./(2*x+b));
xi = @(zeta) zeta.^2 +b.*zeta + c;
epsi_c = A.*(log(x.^2 ./ xi(x)) + 2*b.*eta./Q - (b*x0./xi(x0)) .* (log((x-x0).^2 ./ xi(x)) + (2.*(2*x0+b).*eta)/Q));

Vc = epsi_c - A/3 * (c*(x-x0)-b*x*x0)/(xi(x).*(x-x0));
Vxc_rho = Vx + Vc';

for mu=1:M
    for nu=1:M
        Vxc(mu,nu) = sum(basisfxn(:,mu).*basisfxn(:,nu).*Vxc_rho);
    end
end
Exc =sum((epsi_x + epsi_c').*rho);

end

