%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
% sample driver script for CMOD5.m
%
clear

% compute a single sigma^0 value for a particular geometry and wind vector
spd=10;  % wind speed
dir=30;  % wind direction
az=45;   % radar azimuth angle
inc=30;  % incidence angle

sig=CMOD5([spd, dir, az, inc])

% compute several sigma^0 values at once

spds=[10, 10, 10];
dirs=[30, 30, 30];
azs=[0, 45, 90];
incs=[30, 30, 30];

sigs=CMOD5([spds', dirs', azs', incs'])'
