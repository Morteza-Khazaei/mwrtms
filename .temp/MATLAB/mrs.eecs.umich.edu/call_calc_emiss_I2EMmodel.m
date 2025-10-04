%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

sig1 = 0.88/100; L1 = 10.0/100;
sig2 = 3.00/100; L2 = 20.0/100;
sp = 1; %exponential correl func

theta_d = 40.0; %incidence angle

fr = 1.5; %freq
k = 2*pi*fr/30;
ks1 = k*sig1;
ks2 = k*sig2;
kl1 = k*L1;
kl2 = k*L2;

mv = 0.001: 0.01: 0.35;

np = length(mv);

temp  = 23; % temperature in C

rho_b = 1.7; % g/cm3
Sand = 30.6 / 100;  
Clay = 13.5 / 100;


ev = zeros(np,2);
eh = zeros(np,2);

for n = 1: np
    mv(n)
  [epsr epsi] = RelDielConst_Soil(fr,temp, rho_b, mv(n), Sand, Clay);
  
  er = epsr -1i*epsi;  % permittivity of soil - complex

  [ev(n, 1) eh(n,1)] = Calc_emissivity_I2EMmodel(fr, sig1, L1, theta_d, er, sp);
  [ev(n, 2) eh(n,2)] = Calc_emissivity_I2EMmodel(fr, sig2, L2, theta_d, er, sp);
 
end


figure(3)
plot(mv, ev(:,1),'-r', mv, ev(:,2), ':b')
xlabel('m_v')
ylabel('e_v')
legend('smooth', 'rough')
%axis([0 0.35 0 1])
grid
    
figure(4)
plot(mv, eh(:,1),'-r', mv, eh(:,2), ':b')
xlabel('m_v')
ylabel('e_h')
legend('smooth', 'rough')
%axis([0 0.35 0 1])
grid







