%--------------------------------------------------------------------------
%
% Integrator: Compute the integration with the selected model and on a
%             specific time range
%
% Inputs:
%   start_time      Initial time of the integration (in datetime format ie : 'dd-MM-yyyy')
%   end_time        Last time of the integration (in datetime format ie : 'dd-MM-yyyy')
%   time_step       timestep for the integration in sec 
%   x0              Initial state variables of the S/C (6*1 array) in [m] and [m/s] - ICRF
%                   centered on the Moon (for moon satellite)
%   all the rest    Necessary for the other functions to be used
% 
% Output:
%   tResult    		Integration proper time
%   xResult         State variable of the satellite after integration [m]
%                   for position and [m/s] for velocity 
%
% Last modified:   17/aug/2023   Louis Carton
%
%--------------------------------------------------------------------------


function [tResult, xResult] = Integrator_parfor(x0,options,real_start_time,C,S,C_EGM2008,S_EGM2008,radius_moon,mu_moon,N,AuxParam,SatParam,delta_t)

%MATLAB
%60sec of delta-t = 1440
%10sec of delta-t = 8640
%1sec of delta-t = 86400

time_step = 86400/delta_t;

tspan = linspace(0,86400*AuxParam.duration_prop,(time_step*AuxParam.duration_prop)+1); 

%tspan = [0,86400*AuxParam.duration_prop]; %do not care about the delta-t -> let Matlab decide

[tResult,xResult] = ode78(@(t,x_sc) Accel_Total_parfor(t,x_sc,real_start_time,C,S,C_EGM2008,S_EGM2008,radius_moon,mu_moon,N,AuxParam,SatParam),tspan,x0,options);
%When I am going to need the whole history of the S/C that is contained in
%xResult
end


