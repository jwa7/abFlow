% boysF   Boys function
%
%  y = boysF(m,T)
%
% boysF evaluates the Boys function of order m at argument T. The Boys function
% of order m is defined as
%
%     F_m(T) = int_0^1 t^(2m) exp(-T t^2) dt
%
% boysF can calculate F_m(T) for an array of m or for an array of T.
%
% Examples:
%
% >> boysF(2,0.3)  % scalar m and scalar T
% ans =
%     0.1618
% >> boysF([0 1 2],0.3) % array of m, scalar T
% ans =
%     0.9084    0.2793    0.1618
% >> boysF(2,0.1:0.1:0.5) % scalar m, array of T
% ans =
%     0.1863    0.1735    0.1618    0.1509    0.1408

function y = boysF(m,T)

% Evaluate in terms of the lower incomplete gamma function.
% Warning: Compared to the common definition of the lower
% incomplete gamma function (Wikipedia, Digital Library of
% Mathematical Functions, Wolfram Functions), Matlab's
% gammainc() includes a normalization/regularization factor,
% and the order of arguments is reversed.
mp = m+1/2;
y = (gammainc(T,mp,'lower').*gamma(mp))./(2*T.^mp);

% Limit for T -> 0
threshold = 1e-13;
if numel(T)>1
  idx = abs(T)<threshold;
  y(idx) = 1/(2*m+1);
else
  if abs(T)<threshold
    y = 1./(2*m+1);
  end
end
