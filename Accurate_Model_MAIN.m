% Step by step 

%Set up the constant + access to files of interest
clearvars; close all; 
clear all
clc 

%------------------------------------------------------------------------
%----------- Mac ------------------
userdir ='/Users/louis/Documents/MATLAB/';
%  Add specified folder to the top of the search path
addpath([ userdir , '/mice_mac/lib' ])
addpath([ userdir , '/mice_mac/src/mice' ])
% Construct a meta kernel, “standard.tm”, which will be used to load the needed
% generic kernels: “naif0011.tls”, “de421.bsp”, and “pck00010.tpc”.
% Load the generic kernels using the meta kernel, and a Cassini spk.
cspice_furnsh([ userdir , '/spice_data/de421.bsp' ]);                       % Planetary ephemerides
cspice_furnsh([ userdir , '/spice_data/pck00010.tpc' ]);                    % Planet orientation and radii
cspice_furnsh([ userdir , '/spice_data/naif0011.tls' ]);
cspice_furnsh([ userdir , '/spice_data/gm_de431.tpc' ]);
cspice_furnsh([ userdir , '/spice_data/moon_080317.tf' ]);
cspice_furnsh([ userdir , '/spice_data/moon_pa_de421_1900-2050.bpc' ]);

%---------------------------------------------------------

%------------------------------------------------------------------------
% On Windows
% userdir = getenv("C:/Program Files/");
% 
% % Load SPICE kernel files into Matlab
% kernels = {'spice_data/naif0012.tls','spice_data/pck00010.tpc','spice_data/gm_de431.tpc','spice_data/moon_080317.tf','spice_data/moon_pa_de421_1900-2050.bpc','spice_data/de430.bsp'};
% cspice_furnsh(kernels);
%------------------------------------------------------------------------

%Load the necessary data
%global C S radius_moon mu_moon N AuxParam SatParam
%Get the data that are constant throughout the integration
load mu_moon.mat; load radius_moon.mat; load S_grgm600a.mat; load C_grgm600a.mat; load Normalization_factors.mat;
load C_EGM2008.mat; load S_EGM2008.mat; 


%Select what forces you want to include in the integration + duration of propagation: 
AuxParam = struct('degree',10,'order',10, ...
    'Sun',0,'Earth',1,'degree_Earth',0,'order_Earth',0,'Planets',0,...
    'SRP',1,...
    'Albedo',0,...
    'Solid_Tides',0,...
    'Relativity',0);
AuxParam.duration_prop = 100;%DAYS %duration of the propagation in days!!

SatParam = struct('Area',1,'Mass',1018,'Cr',1.8);
%Area = 12m^2 from the official NASA Layout | mass = 1830; %kg Wet mass | Cr = 1.5; %see "Satellite Orbits" p78 for more info 
    
    
% Data from LRO the 10/03/2012 in the ICRF Moon-centered frame (MI) - Initial conditions + options of integration
r0 = [1.849680295540107E+03 -5.461316594356818E+01 4.501462442022810E+01]; %km
v0 = [5.308416235529445E-03 -6.360739643932042E-01 1.498335232373685];  %km/s
x0 = [r0 v0]'*1e3; %in m and m/s

% Data from LRO the 30/09/2021 in the ICRF Moon-centered frame (MI) 
% r0 = [-1.622249719003224E+03 -8.395251489285182E+02 2.139813531213992E+02]; %km
% v0 = [4.794420452920704E-01 -5.558308717823630E-01 1.455589707948305];  %km/s
% x0 = [r0 v0]'*1e3; %in m and m/s


%Options of Integration
options = odeset('RelTol',1e-12,'AbsTol',1e-14);
start_date_UTC = juliandate(datetime('2012-03-10 00:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS'));
sat = 'LRO';


%% Integration
clc
tic
%[tResult, xResult] = Integrator(real_tspan,fake_tspan,x0,options);
[tResult, xResult] = Integrator_parfor(x0,options,start_date_UTC,C,S,C_EGM2008,S_EGM2008,radius_moon,mu_moon,N,AuxParam,SatParam,2);
toc
%tResult_Gregorian = datetime(tResult,'convertfrom','juliandate');
text_save = ['Accurate_Model_',sat,'_RESULTS_',num2str(AuxParam.degree),'x',num2str(AuxParam.order),'_GSH_Earth_',num2str(AuxParam.degree_Earth),'x',num2str(AuxParam.order_Earth)];
save('/Users/louis/MATLAB-Drive/RESEARCH/Accurate_Model_MATLAB_Codes/Results - Saved Data/text_save',"xResult","tResult")

%% INTEGRATOR DATA (if applicable)


load Accurate_Model_LRO_Results_10x10_Earth_5x5_5d_10_Mar_2012_V2

pos_Integrator = sqrt(xResult(:,1).^2 + xResult(:,2).^2 + xResult(:,3).^2);

% load Sensitivity_Model_Reduction_LRO_RESULTS.mat
% pos_J2 = sqrt(XRESULTS(1,:,2).^2 + XRESULTS(2,:,2).^2 + XRESULTS(3,:,2).^2);
% % pos_5x5 = sqrt(XRESULTS(1,:,7).^2 + XRESULTS(2,:,7).^2 + XRESULTS(3,:,7).^2);
% error_pos = (pos_Integrator - pos_10x10 ) ;
% close(figure(1));figure(1)
% plot(tResult,error_pos)
% xlabel('Time - UTC')
% ylabel('Total Distance Error - [m]')
%% HORIZON DATA
%
%Import the data of LRO from 10/Mar/2012 - 11/Mar/2012 in MI Frame - 
% file_HORIZON = 'C:\Users\louis\MATLAB Drive\RESEARCH\JPL_Horizons_MATLAB_Codes\Horizons_Data-Sets\UTC Dates\-85_500@301_FRAME_VECTORS_2012-Mar-10 UTC_2012-Mar-11';
% SAT_HORIZON = Open_horizons_results_modified(file_HORIZON);


%% GMAT DATA

%Import the data of LRO from 10/Mar/2012 - 11/Mar/2012 in MI Frame - 
file_GMAT = 'LRO_MoonICRF_UTC_2012-03-10_2012-03-20_data_LP165P_2BP';

data = readlines(file_GMAT);
data(end) = []; %remove last line that is empty (for whatever reason)
data = split(data,' ');
data_1 = data(:,1:6);

data_2 = append(data(:,9),'-',data(:,8),'-',data(:,7),' ',data(:,10));
date_GREGORIAN = datetime(data_2,'Format','yyyy-MMM-dd HH:mm:ss.SSS');
SAT_GMAT.cart = str2double(data_1);
SAT_GMAT.dates_GREGORIAN = date_GREGORIAN;


%% Error in position 
% SAT_GMAT.cart = SAT_GMAT.cart(1:172801,:);
% SAT_GMAT.dates_GREGORIAN = SAT_GMAT.dates_GREGORIAN(1:172801);
pos_Integrator = sqrt(xResult(:,1).^2 + xResult(:,2).^2 + xResult(:,3).^2);
pos_GMAT = sqrt(SAT_GMAT.cart(:,1).^2 + SAT_GMAT.cart(:,2).^2 + SAT_GMAT.cart(:,3).^2)*1e3;
error_pos = abs(pos_Integrator - pos_GMAT) ;
close(figure(1));figure(1)
subplot(2,1,1)
hold on
plot(SAT_GMAT.dates_GREGORIAN,pos_Integrator)
plot(SAT_GMAT.dates_GREGORIAN,pos_GMAT)
legend('Integrator','GMAT')
xlabel('Time - UTC')
ylabel('Distance - [m]')
subplot(2,1,2)
plot(SAT_GMAT.dates_GREGORIAN,error_pos)
xlabel('Time - UTC')
ylabel('Total Distance Error - [m]')

%% Positon error component by component Integrator vs GMAT
%Stop comparing with Horizon now that we have a very accurate model on GMAT
close(figure(2));figure(2)
subplot(3,1,1)
hold on 
% plot(tResult_Gregorian,xResult(:,1)*1e-3)
plot(SAT_GMAT.dates_GREGORIAN,xResult(:,1)*1e-3-SAT_GMAT.cart(:,1))
%plot(SAT_GMAT.dates_GREGORIAN,SAT_GMAT.cart(:,1))
legend('Integrator','GMAT')
xlabel('Time - UTC')
ylabel('X - ICRF [km]')
subplot(3,1,2)
hold on 
%plot(tResult_Gregorian,xResult(:,2)*1e-3)
plot(SAT_GMAT.dates_GREGORIAN,xResult(:,2)*1e-3-SAT_GMAT.cart(:,2))
%plot(SAT_GMAT.dates_GREGORIAN,SAT_GMAT.cart(:,2))
legend('Integrator','GMAT')
xlabel('Time - UTC')
ylabel('Y - ICRF [km]')
subplot(3,1,3)
hold on 
%plot(tResult_Gregorian,xResult(:,3)*1e-3)
plot(SAT_GMAT.dates_GREGORIAN,xResult(:,3)*1e-3-SAT_GMAT.cart(:,3))
%plot(SAT_GMAT.dates_GREGORIAN,SAT_GMAT.cart(:,3))
legend('Integrator','GMAT')
xlabel('Time - UTC')
ylabel('Z - ICRF [km]')

%% Pos + Vel of the Integration alone
close(figure(3));figure(3)
subplot(6,1,1)
hold on 
plot(tResult,xResult(:,1)*1e-3)
xlabel('Time - [s]')
ylabel('X - ICRF [km]')
subplot(6,1,2)
hold on 
plot(tResult,xResult(:,2)*1e-3)
xlabel('Time - [s]')
ylabel('Y - ICRF [km]')
subplot(6,1,3)
hold on 
plot(tResult,xResult(:,3)*1e-3)
xlabel('Time - [s]')
ylabel('Z - ICRF [km]')
subplot(6,1,4)
hold on 
plot(tResult,xResult(:,4)*1e-3)
xlabel('Time - [s]')
ylabel('VX - ICRF [km/s]')
subplot(6,1,5)
hold on 
plot(tResult,xResult(:,5)*1e-3)
xlabel('Time - [s]')
ylabel('VY - ICRF [km/s]')
subplot(6,1,6)
hold on 
plot(tResult,xResult(:,6)*1e-3)
xlabel('Time - [s]')
ylabel('VZ - ICRF [km/s]')

%%
close(figure(4));figure(4)
plot3(SAT_GMAT.cart(:,1),SAT_GMAT.cart(:,2),SAT_GMAT.cart(:,3))
