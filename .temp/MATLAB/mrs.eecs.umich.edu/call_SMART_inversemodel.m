%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

%-- start with the forward model
eps = 8.00 - 1i * 1.53; % corresponding to mv = 0.16

f = 7; % frequency in GHz
s = 0.37; % rms height (cm)

theta = 30:1:70;
nt = length(theta);

for ii = 1:nt
[sig_0_vv(ii) sig_0_hh(ii)] = SMART_ForwardModel(eps,theta(ii),s,f);
end

figure(1)
plot(theta, sig_0_vv, '-r', theta, sig_0_hh, '-b')
grid
xlabel('\theta (deg)')
ylabel('\sigma^o  (dB)')
legend('vv', 'hh')

axis([ 30 70 -50 20])

%--- now the inverse model

for ii = 1:nt
[eps_pr(ii) ss(ii)] = SMART_InverseModel(sig_0_vv(ii), sig_0_hh(ii),theta(ii),f);
end

figure(2)
subplot(2,1,1)
plot(theta, ss, '-r') 
xlabel('\theta (deg)')
ylabel('s (cm)')
grid

subplot(2,1,2)
plot(theta, eps_pr, '-b')
xlabel('\theta (deg)')
ylabel('\epsilon ')
grid
