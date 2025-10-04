%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

%-- To Generate plots for complex permittivity at two different
%temperatures

f = 0.1:0.1:100; %frequency in GHz (input)
t1= -1 ; % temperature in degrees C

[epsr1, epsi1]= RelDielConst_PureIce(t1,f);


t2 = -15;
[epsr2, epsi2]= RelDielConst_PureIce(t2,f);

loglog(f, epsi1, f, epsi2)
grid

