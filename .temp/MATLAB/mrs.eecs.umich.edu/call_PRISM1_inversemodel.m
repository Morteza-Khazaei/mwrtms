%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

%--start by generating the data to use in the inversion part
f = 9.5; % frequency in GHz
s = 0.0032; % rms height (m)

mv0 = 0.15; %volumetric moisture
S = 0.306; C = 0.135; rho_b = 1.7; %loam soil parameters
[er ei] = RelDielConst_Soil(f,23, rho_b, mv0, S, C); %get the dielectric
eps2 = er -1i*ei;
abs(eps2)

theta = 5:1:80;
nt = length(theta);

sig_0_vv = zeros(nt,1); %intialize arrays
sig_0_hv = sig_0_vv; sig_0_hh = sig_0_hv;

for ii = 1:nt
[sig_0_vv(ii) sig_0_hh(ii) sig_0_hv(ii)] = PRISM1_ForwardModel(eps2,theta(ii),s,f);
end

figure(1)
plot(theta, sig_0_vv, '-r',  theta, sig_0_hh, '-b', theta, sig_0_hv, '-g')
xlabel('incidence angle')
ylabel('\sigma^o (dB)')
legend('vv-pol', 'hh-pol', 'hv-pol',3)
grid

%-----Use the data in the inversion step ---------------------------------
eps = zeros(nt, 1); s = eps; mv = eps; %intialize arrays

for n = 1: nt    
    [eps(n) s(n) mv(n)] = PRISM1_InverseModel(sig_0_vv(n),sig_0_hh(n), sig_0_hv(n), theta(n), f);
end

s = s * 100; % transform to cm scale

figure(2)
subplot(3,1,1)
plot(theta, eps)
xlabel('\theta (deg)')
ylabel('|\epsilon|')
grid

subplot(3,1,2)
plot(theta, s)
xlabel('\theta (deg)')
ylabel('|s (cm)')
axis([ 0 80 0 max(s)*1.1])
grid

subplot(3,1,3)
plot(theta, mv)
xlabel('\theta (deg)')
ylabel('m_v (g/cm^3)')
axis([ 0 80 0 0.5])
grid
