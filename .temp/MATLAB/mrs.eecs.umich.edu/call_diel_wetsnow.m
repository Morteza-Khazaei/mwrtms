%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

%-- calculate permittivity of wet snow for four different moisture contents
f = 0:0.1:40; %frequency in GHz (input)
mv1 = 12; % volumetric moisture content (percentage)
Ps = 0.24; %snow density (g/cm3)

[epsr1, epsi1] = RelDielConst_WetSnow(Ps, mv1, f );

mv2 = 8; 
[epsr2, epsi2] = RelDielConst_WetSnow(Ps, mv2, f );

mv3 = 4; 
[epsr3, epsi3] = RelDielConst_WetSnow(Ps, mv3, f );
mv4 = 0; 
[epsr4, epsi4] = RelDielConst_WetSnow(Ps, mv4, f );

% subplot(2,1,1)
plot(f, epsr1, f, epsr2, f, epsr3, f, epsr4)
legend('mv=12','8', '4', '0', 1)
xlabel('frequency (GHz)')
ylabel('\epsilon_r')
grid

% subplot(2,1,2)
% plot(f, epsi1, f, epsi2, f, epsi3, f, epsi4)
% xlabel('frequency (GHz)')
% grid

