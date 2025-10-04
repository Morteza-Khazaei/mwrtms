%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear

% select a true wind vector
true_speed=10.0; % true wind speed in m/s
true_dir=120;    % true wind dir in deg

% define radar geometry
azimuth=[45,90,135]';  % beam azimuth angles
incidence=[30,30,30]'; % beam incidence angles

% compute the sigma-0 values corresponding to measurement geometry and wind
z=zeros(size(azimuth));
arg=[true_speed+z,true_dir+z,azimuth,incidence];
sigs=CMOD5(arg)';

% define a particular direction to search
try_dir=25;

% define global variable to pass geometry
global search_geom
% set global variable to pass function arguments                            
search_geom=[try_dir+z,azimuth,incidence,sigs];
 
% try a particular speed
spd_trial=7.0;

% compute objective function value (CMOD5)
LLS=dirmin(spd_trial)

% For a fixed direction, find the wind speed that minimizes LLS
% set option for matlab optimization routine to speed search
options=optimset('TolX',0.1); % only need 0.1 m/s speed precision

% call matlab optimization routine fminbnd
[speed_opt final_LLS]=fminbnd(@dirmin,1,50,options)
