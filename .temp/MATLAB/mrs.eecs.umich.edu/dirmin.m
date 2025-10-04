%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
% Code: 16.3: Utility for computing L_LS objective function 

% Description:  For a particular wind direction, radar geometry, and 
% sigma^0 values, the code evaluates L_LS using CMOD5.m Code is designed 
% to be called by a Matlab minimization routine to support a ravine-search 
% for wind retrieval.
% 
% Input Variables: (vectors, one entry for each beam)
% U: wind speed to try
% phi_wind: wind direction  to try
% phi_radar: radar illumination azimuth angle relative to north (deg)
% theta: incidence angle [15 <= theta <= 60] (deg)
% sigma^0: observed C-band VV-pol backscatter value
% 
% Output Products:
% L_LS:  Least-squares objective function
% 
% Book Reference: Chapter 16 

function LS=dirmin(spd)
%
% function to be used with fminbnd to compute the least-squares 
% objective function at a particular speed using CMOD5 and geometry 
% stored in the global variable search_geom

global search_geom

% generate model function input argument for fixed geometry
arg=[spd+zeros(size(search_geom(:,1))) search_geom(:,1:3)];

% call model function
sigs=CMOD5(arg);

% compute least-squares error
LS=sum((sigs-search_geom(:,4)').^2);
end

