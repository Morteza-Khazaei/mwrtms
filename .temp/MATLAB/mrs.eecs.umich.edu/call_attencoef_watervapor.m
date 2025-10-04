%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

%-- sample driver code (simple case):

t = -10:0.5: 50; % temperature in degree C.
P = 1013; % pressure in mbar
f = 137.8; % frequency in GHz
rho_0 = 10 ; % water vapor density (g/cm3)

np = length(t);
kappaH2O = zeros(np,1);

for ii = 1: np
    [kappaH2O(ii)]=AttenCoef_WaterVapor(t(ii),P,rho_0,f);
end

plot(t, kappaH2O)
grid



% 
% %-- The following code reproduces the Fig. 2a in the paper by Liebe,
% %Hufford, and Cotton, titled "Propagation Modeling of Moist Air and
% %Suspended Water/ICE particles at Frequencies Below 1000 GHz". 
% 
% P = 1013; % Pressure in mbar
% 
% t = 7: 2: 50; % temperature in degrees C.
% 
% f = 137.8; % freq in GHz
% u = 80;
% theta = 300 ./ (t + 273);
% es = 2.408e11 .* theta.^5 .* exp(-22.644 .*theta);
% e = u * es ./100;
% 
% rho_0 = 0.7223 .* e .* theta;
% 
% np = length(t);
% kappaH2O = zeros(np,1);
% 
% 
% for ii = 1: np
%     [kappaH2O(ii)]=AttenCoef_WaterVapor(t(ii),P,rho_0(ii),f);
% end
% 
% plot(t+273, kappaH2O/f/0.182)
% grid
