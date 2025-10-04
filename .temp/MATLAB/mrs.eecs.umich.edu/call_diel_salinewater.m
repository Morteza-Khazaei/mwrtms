%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

%-- To Generate plots for complex permittivity of saline water.

f = 0:1:1000; %frequency in GHz
t= 20; % temperature in degree C
S = 32.54; % percentage salinity


[epsr, epsi]= RelDielConst_SalineWater(t,f, S);



subplot(2,1,1)
loglog(f, epsr)
grid

subplot(2,1,2)
loglog(f,epsi)
grid