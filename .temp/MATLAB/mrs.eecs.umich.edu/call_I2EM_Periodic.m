%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%

% This is a computationally demanding program, so it may take several
% minutes for it to execute. 

clear;

%-- list of input parameters

f = 1.9; %1.9; % freq in GHz
eps = 14 -1i*7; %14 - 1i*7; % dielectric constant

%-periodic surface: f(x,y) = -A cos(2 pi y/Gmm)
A = 0.25; %0.25; %Amplitude value in meters
Gmm = 1.0; %1.0; %period in meters

sp = 1; % 1 for exponential and 2 for gaussian correlation functions
s = 0.01; %0.01; % rms height in m
l = 0.041; %0.041; % correl length (m)


%--------------------------------------------
%-- vary theta but fix phi
ph_deg = 88.5; %88.5; %phi = 0 is when radar is parallel to rows
                        %phi = 90 is when radar is perpendicular to rows
n_ph = length(ph_deg);

n_th = 71;  % number of angles
th_start = 1; %start angle in deg
th_stop =71; % stop angle in deg

the_deg = linspace(th_start,th_stop, n_th); %create the vector of incident angles


sig_o_vv = zeros(n_ph, n_th); % initialize arrays
sig_o_hh = sig_o_vv;
sig_o_vh = sig_o_vv;

[sig_o_vv, sig_o_hh, sig_o_vh] = I2EM_Periodic(the_deg,ph_deg ...
    , eps, f, s, l, sp, Gmm, A);

    %To maintain reasonable execution time for this code, calculation of
    %sigma_0_hv is switched off here. However, it can be activated with
    %ease (check the function called Calc_scatt_coeff listed below).

sig_o_vv_dB= 10*log10(sig_o_vv); % case into dB scale
sig_o_hh_dB= 10*log10(sig_o_hh);
sig_o_vh_dB= 10*log10(sig_o_vh);

figure(1) 
plot(the_deg, sig_o_vv_dB(1,:), the_deg, sig_o_hh_dB(1,:) )
set(gca, 'FontSize', 12, 'FontWeight', 'b')
xlabel('Incidence Angle \theta (deg)','FontSize',14, 'FontWeight','b')
ylabel('\sigma^o (dB)', 'FontSize',14, 'FontWeight','b')
%title('Volume Scattering - Lossy Sand Particles', 'FontSize',16, 'FontWeight','b')
ch = get(gca, 'Children');
set(ch(2),'LineWidth', 1.5, 'LineStyle','-','Color','r')
set(ch(1),'LineWidth', 1.5, 'LineStyle','-','Color','b')
%axis([0 50 -40 5])
legend ('vv', 'hh', 3)
grid






% %--------------------------------------------
% %-- vary phi but fix theta
% ph_deg = 0:2:90; %phi = 0 is when radar is parallel to rows
%                         %phi = 90 is when radar is perpendicular to rows
% n_ph = length(ph_deg);
% 
% the_deg = 60;
% n_th = 1;
% 
% sig_o_vv = zeros(n_ph, n_th); % initialize arrays
% sig_o_hh = sig_o_vv;
% sig_o_vh = sig_o_vv;
% 
% for jj = 1: n_ph % loop through all angle combinations of interest to generate the data
%     for ii = 1: n_th
%         the = the_deg(ii) ;
%         ph = ph_deg(jj)
%         [sig_o_vv(jj,ii), sig_o_hh(jj,ii), sig_o_vh(jj,ii)] = I2EM_Periodic(the,ph ...
%          , eps, f, s, l, sp, Gmm, A);
% 
%     %To maintain reasonable execution time for this code, calculation of
%     %sigma_0_hv is switched off here. However, it can be activated with
%     %ease (check the function called Calc_scatt_coeff listed below).
% 
%     end
% end
% 
% sig_o_vv_dB= 10*log10(sig_o_vv); % case into dB scale
% sig_o_hh_dB= 10*log10(sig_o_hh);
% sig_o_vh_dB= 10*log10(sig_o_vh);
% 
% figure(2) 
% plot(ph_deg, sig_o_vv_dB(:,1), ph_deg, sig_o_hh_dB(:,1) )
% set(gca, 'FontSize', 12, 'FontWeight', 'b')
% xlabel('Azimuth Angle \theta (deg)','FontSize',14, 'FontWeight','b')
% ylabel('\sigma^o (dB)', 'FontSize',14, 'FontWeight','b')
% %title('Volume Scattering - Lossy Sand Particles', 'FontSize',16, 'FontWeight','b')
% ch = get(gca, 'Children');
% set(ch(2),'LineWidth', 1.5, 'LineStyle','-','Color','r')
% set(ch(1),'LineWidth', 1.5, 'LineStyle','-','Color','b')
% %axis([0 50 -40 5])
% legend ('vv', 'hh', 3)
% grid




