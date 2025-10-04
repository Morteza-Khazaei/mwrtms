%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

%-- sample driver code (simple case): reproduces fig 8-7

t = 20; % temperature in degree C.
P = 1013; % pressure in mbar
f = 1:1:300; % frequency in GHz
rho_0 = 7.5 ; % water vapor density (g/cm3)

np = length(f);
kappag = zeros(np,1);

for ii = 1: np
    [kappag(ii)]=AttenCoef_Gaseous(t,P,rho_0,f(ii));
end

figure(1)
loglog(f, kappag)
grid




% %-- The following code reproduces the Fig. 2a in the paper by Liebe,
% %Hufford, and Cotton, titled "Propagation Modeling of Moist Air and
% %Suspended Water/ICE particles at Frequencies Below 1000 GHz". 
% 
% P = 1013; % Pressure in mbar
% 
% t = 20; % temperature in degree C.
% f = 1:1:1000; % frequency in GHz
% u = 0;
% theta = 300 ./ (t + 273);
% es = 2.408e11 .* theta.^5 .* exp(-22.644 .*theta);
% e = u * es ./100;
% 
% rho_0 = 0.7223 .* e .* theta;
% 
% np = length(f);
% kappaO2 = zeros(np,1);
% 
% 
% for ii = 1: np
%     [kappaO2(ii)]=AttenCoef_Oxygen(t,P,rho_0,f(ii));
% end
% 
% semilogy(f, kappaO2 )
% grid
