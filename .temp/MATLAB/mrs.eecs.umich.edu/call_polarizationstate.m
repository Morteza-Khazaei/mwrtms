%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

a_y = 1;

a_x = 0.01:0.01:100;

delta = -40;

np = length(a_x);
psi = zeros(np,1);
chi = zeros(np,1);

for ii = 1: np
[psi(ii) chi(ii)] = PolarizationState(a_x(ii), a_y, delta);
end

subplot(2,1,1)
semilogx(a_x, psi)
grid

subplot(2,1,2)
semilogx(a_x, chi)
grid