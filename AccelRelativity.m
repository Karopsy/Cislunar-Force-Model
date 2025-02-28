%--------------------------------------------------------------------------
%
% AccelRelativity: Computes the perturbational acceleration due to relativistic
%              effects (ref : "General relativistic effects acting on the
%              orbit of Galileo satellites" - see bibliography)
%
% Inputs:
%   r_sat           Satellite position vector (ICRF) [m]
%   v_sat           Satellite velocity vector (ICRF) [m/s]
%   mu              Gravitational coefficient of primary body [m^3/s^2]
% 
% Output:
%   a_REL    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   9/aug/2023   Louis Carton
%
%--------------------------------------------------------------------------

function a_REL = AccelRelativity(r_sat,v_sat,mu)

c = 299792458; % speed of light in m/s
r = norm(r_sat);
v = norm(v_sat);
beta = 1; %relativistic coefficient
gamma = 1; %relativistic coefficient

a_REL = abs((mu/(c^2*r^3)).*((2*(beta+gamma)*mu/r - gamma*v^2).*r_sat + 2*(1+gamma)*v^2.*v_sat)); 
    %+ (1+gamma)*(mu/(c^2*r^3))*((3/2*r^2)*cross(r_sat,v_sat)*dot(r_sat,J)+cross(v_sat,J));

%Line 1 : Schwarzschild correction (~10e-10 m/s^2)
%Line 2 : Lense-Thirring correction (~10e-11 m/s^2) - J is the angular
%momentum per unit mass (i.e cross(r,v)) BUT for the Moon!
end










