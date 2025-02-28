%--------------------------------------------------------------------------
% 
% Conical
%
% Purpose:
%   Computes the fractional illumination of a spacecraft in the 
%   vicinity of the Moon assuming a conical shadow model
% 
% Inputs:
%   r           Spacecraft position vector (ICRF) [m]
%   r_occbody   occulting body (Earth or Moon) (ICRF - rotating frame) [m]
%   r_Sun       Sun position vector (ICRF) [m]
%   
% Output:
%   nu        Illumination factor:
%             nu=0   Spacecraft in Earth shadow
%             0<nu<1 Spacecraft in Penumbra
%             nu=1   Spacecraft fully illuminated by the Sun
%
% Reference:
% Montenbruck O., and Gill E., "Satellite Orbits: Models, Methods, and 
% Applications," Springer Verlag, Heidelberg, Corrected 3rd Printing (2005).
% 
% Last modified:   11/aug/2023   Louis Carton
% 
%--------------------------------------------------------------------------

function nu = Conical(r_sat,r_Sun,r_occbody,occ_body)

R_earth = 6378e3;%m
R_moon = 1738e3;%m
R_sun = 696342e3;%m

if strcmp(occ_body,'Earth')
    R_occbody = R_earth;
elseif strcmp(occ_body,'Moon')
    R_occbody = R_moon;
end

%s_sun = r_Sun - r_occbody;
s = r_sat - r_occbody;
s_mag = norm(s);

d = r_Sun - r_sat;
d_mag = norm(d);

a = asin(R_sun/d_mag);          % eq. 3.85
b = asin(R_occbody/s_mag);        % eq. 3.86
c = acos(-1.0*dot(s,d)/(s_mag*d_mag)); % eq. 3.87

if( c >= (a+b) )            % in Sun light
    nu = 1.0;
elseif( c < abs(a-b) )      % in Umbra
    nu =  0.0;
else                        % in Penumbra 
    x = (c^2+a^2-b^2)/(2*c);                  % eq. 3.93
    y = sqrt(a^2-x^2);
    A = a^2*acos(x/a)+b^2*acos((c-x)/b)-c*y;  % eq. 3.92
    nu = 1.0 - A/(pi*a^2);                    % eq. 3.94
end

end


