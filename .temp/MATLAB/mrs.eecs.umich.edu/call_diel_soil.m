%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
% This code is based on two semi-empirical models, one for the range 0.3 to
% 1.4 GHz (Section 4-8.3 in the book) and another for 1.4 to 37 GHz. The
% two models give different values at 1.4 GHz. Hence, plots as a function
% of frequency may exhibit a discontinuity at 1.4 GHz.

clear;

f = 0.3:.1:18; % frequency in GHz
t= 23; % temperature in degree C

Sand = 50 / 100; % fraction of Sand
Clay = 15 /100; % fraction of Clay
mv = 0.158; % volumetric moisture content
rho_b = 1.7; % soil bulk density g/cm3

nf = length(f);

for ii = 1: nf
[epsr1(ii), epsi1(ii)] = RelDielConst_Soil(f(ii),t, rho_b, mv, Sand, Clay);
end
figure(1)
subplot(2,1,1)
plot(f, epsr1)
xlabel('frequency (GHz)')
ylabel('\epsilon_r')
grid

subplot(2,1,2)
plot(f,epsi1)
xlabel('frequency (GHz)')
ylabel('\epsilon_i')
grid


% f = 6;
% t= 23; % in C
% 
% Sand = 30.6;
% Clay = 13.5;
% mv = 0.0:0.01:0.4;
% rho_b = 1.7; % bulk density g/cm3
% 
% nm = length(mv);
% for ii = 1: nm
% [epsr1(ii), epsi1(ii)] = RelDielConst_Soil(f,t, rho_b, mv(ii), Sand, Clay);
% end
% figure(1)
% subplot(2,1,1)
% plot(mv, epsr1)
% grid
% 
% subplot(2,1,2)
% plot(mv,epsi1)
% grid
% 
% %-----------------------------------------------
% %------------------------------------------------
% 
% clear;
% f = 3.0;
% t = 23;
% Sand = 30.6;
% Clay = 13.5;
% mv = 0.30;
% 
% rho_b = 1.7;
% 
% [epsr1, epsi1] = RelDielConst_Soil(f,t, rho_b, mv, Sand, Clay)


