%--------------------------------------------------------------------------
%
% AccelPointMass: Computes the perturbational acceleration due to a point
%                 mass
%
% Inputs:
%   r_sat          Satellite position vector (ICRF) [m]
%   r_pm           Point mass position vector (ICRF) [m]
%   mu             Gravitational coefficient of point mass [m^3/s^2]
%
% Output:
%   a_PM    	   Acceleration (a=d^2r/dt^2)
%
% Last modified:   9/aug/2023   Louis Carton
%
%--------------------------------------------------------------------------
function a_PM = AccelPointMass(r_sat,r_pm,mu)

% Relative position vector of satellite w.r.t. point mass 
d = r_sat - r_pm;

% Acceleration 
a_PM = -mu * ( d/(norm(d)^3) + r_pm/(norm(r_pm)^3) );

end





