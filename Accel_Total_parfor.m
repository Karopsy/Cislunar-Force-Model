%--------------------------------------------------------------------------
%
% Accel_Total: Computes the acceleration of a Moon orbiting satellite due to 
%    	 - Moon's spherical harmonic gravity field
%        - Moon's gravitationnal centrifugal perturbation
%    	 - gravitational perturbations of the Sun, Earth and Jupiter
%    	 - solar radiation pressure
%	 	 - relativity
%
% Inputs:
%   Mjd_UTC     Modified Julian Date (UTC)
%   X           Satellite state vector in the ICRF coordinate frame
%   Area        Cross-section 
%   mass        Spacecraft mass
%   Cr          Radiation pressure coefficient
%
% Output:
%   dX		    Acceleration (a=d^2r/dt^2) in the ICRF (Inertial) frame
%
% Last modified:   10/Nov/2023   Louis Carton
% 
%--------------------------------------------------------------------------

function dX = Accel_Total_parfor(t,X,real_start_time,C,S,C_EGM2008,S_EGM2008,radius_moon,mu_moon,N,AuxParam,SatParam)

%------------------------------------------------------------------
%--Set-up the variables to compute every acceleration components--- 
%------------------------------------------------------------------

%Get the Time (UTC) at each time step of the integration:
time_UTC = datetime(real_start_time,'convertfrom','juliandate') + t/86400;

%Get the Time (ET) from Time (UTC) (to call SPICE for 3rd body pert)
if isa(time_UTC,'datetime') == 1
    time_UTC = char(time_UTC);
end
time_ET = cspice_str2et(time_UTC);

%Split the vector state
r_MI = X(1:3);
v_MI = X(4:6);

%Compute the rotation matrix from MoonInertial to PrincipalAxis (body-fixed frame) 
rot_mat_M_BF = MI2PA_ET(time_ET); %* PA2ME;

%Compute the vector state in Body Fixed frame (PA or ME)
r_BF = rot_mat_M_BF * r_MI;
%From the vector state in BF frame - compute the spherical components
% r = norm(r_BF);        % distance 
% lat = asin(r_BF(3)/r); %from spherical coordinates
% lon = atan2(r_BF(2),r_BF(1)); %from spherical coordinates

% moon_elevation = Moon_elevation(lat,lon,'rad'); %Integrated in the event
% function @stabilityEvents_2

%Compute the rotation matrix from ICRF (inertial) to Earth Body-Fixed frame
rot_mat_E_BF = ICRF2Earth_BF(time_ET);


degree = AuxParam.degree;
order = AuxParam.order;
Area = SatParam.Area;
Mass = SatParam.Mass;
Cr = SatParam.Cr;
degree_Earth = AuxParam.degree_Earth;
order_Earth = AuxParam.order_Earth;


%------------------------------------------------------------------
%--------------------Compute the Accelerations---------------------
%------------------------------------------------------------------

%Get the positions of Sun and Earth as they are used in various acceleration elements
if (AuxParam.Sun == 1 || AuxParam.Earth == 1 || AuxParam.SRP == 1 || AuxParam.Solid_Tides == 1 || AuxParam.Albedo == 1)
    %SUN
    X_sun = cspice_spkezr('Sun',time_ET,'J2000','NONE','301');
    r_sun = X_sun(1:3)*1e3;
    mu_sun = cspice_bodvrd( 'SUN', 'GM', 1 )*1e9;
    %EARTH
    X_earth = cspice_spkezr('399',time_ET,'J2000','NONE','301');
    r_earth = X_earth(1:3)*1e3;
    mu_earth = cspice_bodvrd( 'EARTH', 'GM', 1 )*1e9;
    r_earth_MoonFixed = rot_mat_M_BF * r_earth;
end

%-----------------------------------------
%Gravity Spherical Harmonics for the Moon:
%-----------------------------------------

%Check if we want to include Solid Tides
if (AuxParam.Solid_Tides == 1)
    %Get the necessary elements:
    %Sun
    r_sun_MoonFixed = rot_mat_M_BF * r_sun;

    %Get the Coefficients correction value
    [C_Corr, S_Corr] = Correction_Solid_Tides(degree,order,r_earth_MoonFixed,mu_earth,r_sun_MoonFixed,mu_sun,mu_moon,radius_moon);

    C_wtides = C(1:degree+1,1:order+1) + C_Corr;
    S_wtides = S(1:degree+1,1:order+1) + S_Corr;
    
    % Gottlieb formulation (faster)
    a_TOT = Accel_Grav_Harmonic_gottlieb(mu_moon, radius_moon, r_MI, C_wtides, S_wtides, degree, order, rot_mat_M_BF);
    % Cunningham formulation
    %     a_GHAR = Accel_Grav_Harmonic_parfor_mex(r_BF,degree,order,C_wtides,S_wtides,mu_moon,radius_moon,N); %Accel in BF frame
    %     a_TOT = rot_mat_M_BF' * a_GHAR; %From BF to Inertial frame
    % end
else
    % if (degree == 0) %2BP
    %     r = norm(r_BF);
    %     ax = -mu_moon*r_BF(1)/r^3;
    %     ay = -mu_moon*r_BF(2)/r^3;
    %     az = -mu_moon*r_BF(3)/r^3;
    %     a_TOT = rot_mat_M_BF' * [ax ay az]';
    % else
    %Gottlieb-Normalized Formulation -- For J2 only need degree=2 and order=0
    a_TOT = Accel_Grav_Harmonic_gottlieb(mu_moon, radius_moon, r_MI, C, S, degree, order, rot_mat_M_BF);
    % end
end




%-----------------------------------------
%Gravity Spherical Harmonics for the Earth:
%-----------------------------------------

if (AuxParam.Earth == 1)

    if (degree_Earth == 0) %Point mass pert
        a_TOT = a_TOT + AccelPointMass(r_MI,r_earth,mu_earth); %Earth pert

    else %Take into account the oblateness of the Earth
        radius_earth = 6378e3; %m
        r_E_ICRF = r_MI - r_earth; 
        a_GHAR_Earth_Inertial = Accel_Grav_Harmonic_gottlieb(mu_earth, radius_earth, r_E_ICRF ,C_EGM2008,S_EGM2008,degree_Earth,order_Earth,rot_mat_E_BF); %Accel in EARTH - Inertial frame !!!!
        a_pert_Earth_Inertial = (-mu_earth/norm(r_earth)^3)*r_earth;
        a_TOT_E_Inertial = a_GHAR_Earth_Inertial + a_pert_Earth_Inertial;
        
        a_TOT = a_TOT + a_TOT_E_Inertial ;
    end
end



%3rd Body Pert:
if (AuxParam.Sun == 1)
    %Sun
    %------ No need to be adapted to parfor--------
    a_TOT = a_TOT + AccelPointMass(r_MI,r_sun,mu_sun); %Sun pert
end

if (AuxParam.Planets == 1)
    % %Mercury
    X_mercury = cspice_spkezr('199',time_ET,'J2000','NONE','301');
    r_mercury = X_mercury(1:3)*1e3;
    mu_mercury = cspice_bodvrd( 'MERCURY', 'GM', 1 )*1e9;
    %Venus
    X_venus = cspice_spkezr('299',time_ET,'J2000','NONE','301');
    r_venus = X_venus(1:3)*1e3;
    mu_venus = cspice_bodvrd( 'VENUS', 'GM', 1 )*1e9;
    %Mars
    X_mars = cspice_spkezr('499',time_ET,'J2000','NONE','301');
    r_mars = X_mars(1:3)*1e3;
    mu_mars = cspice_bodvrd( 'MARS', 'GM', 1 )*1e9;
    %Jupiter
    X_jupiter = cspice_spkezr('599',time_ET,'J2000','NONE','10');
    r_jupiter = X_jupiter(1:3)*1e3;
    mu_jupiter = cspice_bodvrd( 'JUPITER', 'GM', 1 )*1e9;
    %Saturn
    X_saturn = cspice_spkezr('699',time_ET,'J2000','NONE','301');
    r_saturn = X_saturn(1:3)*1e3;
    mu_saturn = cspice_bodvrd( 'SATURN', 'GM', 1 )*1e9;
    %Uranus
    X_uranus = cspice_spkezr('799',time_ET,'J2000','NONE','301');
    r_uranus = X_uranus(1:3)*1e3;
    mu_uranus = cspice_bodvrd( 'URANUS', 'GM', 1 )*1e9;
    %Neptune
    X_neptune = cspice_spkezr('899',time_ET,'J2000','NONE','301');
    r_neptune = X_neptune(1:3)*1e3;
    mu_neptune = cspice_bodvrd( 'NEPTUNE', 'GM', 1 )*1e9;

    a_TOT = a_TOT + AccelPointMass(r_MI,r_mercury,mu_mercury); %Mercury pert
    a_TOT = a_TOT + AccelPointMass(r_MI,r_venus,mu_venus); %Venus pert
    a_TOT = a_TOT + AccelPointMass(r_MI,r_mars,mu_mars); %Mars pert
    a_TOT = a_TOT + AccelPointMass(r_MI,r_jupiter,mu_jupiter); %Jupiter pert
    a_TOT = a_TOT + AccelPointMass(r_MI,r_saturn,mu_saturn); %Saturn pert
    a_TOT = a_TOT + AccelPointMass(r_MI,r_uranus,mu_uranus); %Uranus pert
    a_TOT = a_TOT + AccelPointMass(r_MI,r_neptune,mu_neptune); %Neptune pert
end


%Solar Radiation Pressure
if (AuxParam.SRP == 1) 
    %------ No need to be adapted to parfor--------
    a_TOT = a_TOT + AccelSRP(r_MI,r_earth,r_sun,Area,Mass,Cr,'conical');
end

%Moon Albedo Radiation Pressure (MARP)
if (AuxParam.Albedo == 1)
    a_TOT = a_TOT + AccelAlbedo(r_MI,r_earth,r_sun,Area,Mass,Cr);
end


%General Relativity 
if (AuxParam.Relativity == 1)
    a_TOT =  a_TOT + AccelRelativity(r_MI,v_MI,mu_moon);
end


%------------------------------------------------------------------
%----------------- Final Result to be integrated ------------------
%------------------------------------------------------------------

dX = [v_MI; a_TOT];

end