% Correction 
%--------------------------------------------------------------------------
% Compute the Corrections to bring to the lunar spherical gravity coefficients
%                         -MOON SOLID TIDES-
% 
% Inputs:
% degree        maximum degree
% order         maximum order
% lat           latitude angle (spherical coordinate)
% lon           long angle (spherical coordinate)
% angle_type    'rad' or 'deg'
% r             norm of the position vector in spherical, synodic frame
%
% Outputs:
% Corr_C        Coefficients to add to the GSH coeff (C)
% COrr_S        Coefficients to add to the GSH coeff (S)
%
% Ref: "Analytical Radial Adaptive Method for Spherical Gravity Models" -
% Ahmed Atallah
%
% Last modified:   1/Apr/2024 - Louis Carton
%--------------------------------------------------------------------------

function [Corr_C, Corr_S] = Correction_Solid_Tides(degree,order,r_earth_MoonFixed,mu_earth,r_sun_MoonFixed,mu_sun,mu_moon,radius_moon)


Corr_C = zeros(degree+1,order+1);
Corr_S = zeros(degree+1,order+1);

% Auxiliary quantities - spherical coordinates

%Need the latitude and longitude of the bodies causing the Tides wrt to the
%central body (i.e the moon) 
r_earth = norm(r_earth_MoonFixed);
lat_earth = asin(r_earth_MoonFixed(3)/r_earth); %from spherical coordinates
lon_earth = atan2(r_earth_MoonFixed(2),r_earth_MoonFixed(1)); %from spherical coordinates

r_sun = norm(r_sun_MoonFixed);
lat_sun = asin(r_sun_MoonFixed(3)/r_sun); %from spherical coordinates
lon_sun = atan2(r_sun_MoonFixed(2),r_sun_MoonFixed(1)); %from spherical coordinates

[M_earth,W_earth] = Normalized_Functions_Enhanced_parfor(degree,order,lat_earth,lon_earth,r_earth,radius_moon);
[M_sun,W_sun] = Normalized_Functions_Enhanced_parfor(degree,order,lat_sun,lon_sun,r_sun,radius_moon);

% Love Number for solid tide model: Ref: "Quantification of tidal parameters from Solar System data", 2016
% k(2,0) = 0.02408 ± 0.00045;
% k(2,1) = 0.02414 ± 0.00025;
% k(2,2) = 0.02394 ± 0.00028;
% k(2) = 0.02405; %average
% k(3) = 0.0089;
% k(2,0) = 0.02408;
% k(2,1) = 0.02414;
% k(2,2) = 0.02394;
% k(3,0) = 0.00734;

% Love Number for solid tide model: Ref: GMAT Data 
k(2) = 0.024116;
k(3) = 0;

% Love Number for solid tide model: Ref: "The JPL lunar gravity field to spherical harmonic degree 660 from the GRAIL Primary Mission", 2013

% k(2) = 0.02405;
% k(3) = 0.0089;


%Tides produced by both Earth and Sun
for n=2:3
    for m=0:n
        Corr_C(n+1,m+1) = k(n)/(2*n+1) * ( (mu_earth/mu_moon)*M_earth(n+1,m+1) + (mu_sun/mu_moon)*M_sun(n+1,m+1) );
        Corr_S(n+1,m+1) = k(n)/(2*n+1) * ( (mu_earth/mu_moon)*W_earth(n+1,m+1) + (mu_sun/mu_moon)*W_sun(n+1,m+1) );
    end
end





