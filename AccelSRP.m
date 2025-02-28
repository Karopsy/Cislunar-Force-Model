%--------------------------------------------------------------------------
%
% AccelSRP: Computes the perturbational acceleration due to Solar Radiation
%           Pressure with surface normal to the sun's rays
%                   
% Inputs:
%   r_sat           Satellite position vector (ICRF - body equator of the primary - synodic frame) [m]
%   r_Earth         Earth position vector (ICRF - body equator of the primary - synodic frame) [m]
%   r_Sun           Sun position vector (ICRF - body equator of the primary - synodic frame) [m]
%
%
% Output:
%   a_SRP    		Acceleration (a=d^2r/dt^2)
%
% Ref: "Fundamentals of Astrodynamics and Applications" - David A.Vallado
% 8-43 p585
%
% Last modified:   28/Nov/2023   Louis Carton
%
%--------------------------------------------------------------------------

function a_SRP = AccelSRP(r_sat,r_Earth,r_Sun,Area,mass,Cr,shadow_model)

phi_earth = 1367; %W/m^2 - to be modified if necessary [1360 - 1374]
AU = 1.495978707e11; %m
% In order to have a more accurate value of the solar radiation constant
% that gets to the sat we will use a rule of distance ratio: 

d_satsun = norm(-r_sat + r_Sun);

phi_sat = d_satsun/AU * phi_earth;

P0 = phi_sat/(299792458); %N/m^2 = W*s/m^3 this is the mean solar pressure from the mean solar radiation constant at the S/C distance from the Sun


if contains(shadow_model,'con')
%-----------------------Conical Model---------------------------------------
% For "Accurate_Model" we need to use the shadow function to account for
% the eclipses phenomena for both the Earth and the Moon
    nu_1 = Conical(r_sat,r_Sun,r_Earth,'Earth');
    nu_2 = Conical(r_sat,r_Sun,[0 0 0]','Moon');
    nu = min([nu_1 nu_2]);

elseif contains(shadow_model,'cyl')
%-----------------------Cylindrical Model---------------------------------------
    nu_1 = Cylindrical(r_sat, r_Sun,'Earth');
    nu_2 = Cylindrical(r_sat, r_Sun,'Moon');
    nu = min([nu_1 nu_2]);

elseif contains(shadow_model,'none')
%Exemple for "Force_Evolution_Moon_Environment" we do not use the shadow 
%function as the objective is not to compute an accurate ephem, just to have 
%an idea of the SRP acc value
    nu = 1;
end

% Relative position vector of spacecraft w.r.t. Sun
d = r_sat - r_Sun;

% a_SRP = -nu*Cr*(Area/mass)*P0*d/norm(d)^3*AU^2; % eq. 8-43
a_SRP = nu*Cr*(Area/mass)*P0*d/norm(d); % eq. 8-43


