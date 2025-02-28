%--------------------------------------------------------------------------
% 
% Cylindrical: Computes the fractional illumination of a spacecraft in the 
%              vicinity of the Earth assuming a cylindrical shadow model
% 
% Inputs:
%   r         Spacecraft position vector (ICRF) [m]
%   r_Sun     Sun position vector (ICRF) [m]
%
% Output:
%   nu        Illumination factor:
%             nu=0   Spacecraft is in shadow of the occ body
%             nu=1   Spacecraft is not in the shadow of the occ body
%
% Last modified:   2023/12/5   Louis Carton
%
%--------------------------------------------------------------------------
function nu = Cylindrical(r, r_Sun,occ_body)

R_earth = 6378e3;%m
R_moon = 1738e3;%m

if strcmp(occ_body,'Earth')
    R_occbody = R_earth;
elseif strcmp(occ_body,'Moon')
    R_occbody = R_moon;
end

e_Sun = r_Sun / norm(r_Sun);   % Sun direction unit vector
s     = dot ( r, e_Sun );      % Projection of s/c position 
if ( s>0 || norm(r-s*e_Sun)>R_occbody )
    nu = 1;
else
    nu = 0;
end



