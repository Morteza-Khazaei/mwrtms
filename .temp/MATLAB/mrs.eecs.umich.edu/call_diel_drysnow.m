%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

f = 9.375;%frequency in GHz (input)
t = -18:1:0; % temperature in degree C

Ps = 0.46; % snow density g/cm3
ns = length(t);

for ii = 1: ns 
[epsr(ii), epsi(ii)] = RelDielConst_DrySnow(t(ii), Ps, f);
end

losstangent = epsi ./epsr ;
figure(1)
plot(-t, losstangent)
grid


% f = 0.1:0.1:40;%frequency in GHz (input)
% t = -30; % temperature in degree C
% 
% Ps = 0.45; % snow density g/cm3
% ns = length(f);
% 
% for ii = 1: ns 
% [epsr(ii), epsi(ii)] = RelDielConst_DrySnow(t, Ps, f(ii));
% end
% 
% losstangent = epsi ./epsr;
% 
% figure(1)
%  loglog(f, losstangent)
% % loglog(f, epsi)
% grid

