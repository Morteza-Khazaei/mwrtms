%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

eps = 6.28 - 1i * 1.53; % corresponding to mv = 0.16

f = 9.5; % frequency in GHz
s = 0.00302; % rms height (m)

theta = 5:1:80;
nt = length(theta);

for ii = 1:nt
[sig_0_vv(ii) sig_0_hh(ii) sig_0_hv(ii)] = PRISM1_ForwardModel(eps,theta(ii),s,f);
end


plot(theta, sig_0_vv, '-r',  theta, sig_0_hh, '-b', theta, sig_0_hv, '-g')
xlabel('incidence angle')
ylabel('\sigma^o (dB)')
legend('vv-pol', 'hh-pol', 'hv-pol',3)
grid



