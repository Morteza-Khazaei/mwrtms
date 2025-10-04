%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

eps2 = 3.2 -1i*0.05;
eps3 = 60- 1i*10;

theta_i = 10; 

d = 0:0.01:2; %layer thickness in m

np = length(d); 

f = 1.2; %freq (GHz)


for n = 1:np
[e_v_coh(n) e_h_coh(n)] = CohEmissivity_TwoLayer(eps2, eps3, theta_i, d(n), f) ;

end

figure(1)
plot(d, e_v_coh, '-r', d, e_h_coh,'-b')
xlabel('thickness (m)')
ylabel('Coherent Emissivity')
legend('e_v', 'e_h', 4)
%axis([0 2 0.3 1])
grid