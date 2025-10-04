%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

% this sample code generates the plot in Fig. 6-22 

eps1r = 1;

eps2 = 2.1 -1i * 0.01; % eps of oil film
eps3 = 36 - 1i* 30; % eps of water

f = 20; % freq in GHz

d = 0:0.1: 10; % thickness in mm
d = d / 1000; % in m;
np = length(d);

T0 = 293; % temp in Kalvin

theta = 0; % incidence angle in degrees

%-- brightness temperature of water only (planar boundary)


[eh_w ev_w] = emissivity_PlanarBoundary(eps1r, eps3, theta);
Tb_w = eh_w .*T0;

%-- brightness temp of composite 

for ii = 1: np
    [eh_o(ii) ev_o(ii)] = emissivity_TwoLayerComposite(eps1r, eps2, eps3, d(ii), theta, f);
end

Tb_o = eh_o .* T0;

DeltaTb = Tb_o - Tb_w;

plot(d*1000, DeltaTb)
grid