%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

%-- list of input parameters

f = 1.4; % freq in GHz
eps = 9 - 1i*1.5; % dielectric constant

%-periodic surface: f(x,y) = A cos(2 pi y/Gmm)
A = 0.10; %0.25; %Amplitude value in meters
Gmm = 0.950;%1; %period in meters

sp = 1; % 1 for exponential and 2 for gaussian correlation functions
s = 0.000001; % rms height in m
l = 0.041; % correl length (m)

%-- vary theta but fix phi
ph_deg = 90; %phi = 0 is when radar is parallel to rows
                        %phi = 90 is when radar is perpendicular to rows
n_ph = length(ph_deg);

n_th = 15;  % number of angles
th_start = 0; %start angle in deg
th_stop =70; % stop angle in deg

the_deg = linspace(th_start,th_stop, n_th); %create the vector of incident angles

% the_deg = 10;
% n_th = 1;

ev = zeros(n_ph, n_th); % initialize arrays
eh = ev;

for jj = 1: n_ph % loop through all angle combinations of interest to generate the data
    for ii = 1: n_th
        the = the_deg(ii) 
        ph = ph_deg(jj);
        [ev(jj,ii), eh(jj,ii)] = Emissivity_I2EM_Periodic(the,ph ...
         , eps, f, s, l, sp, Gmm, A);

    end
end

%%

figure(1) 
plot(the_deg, ev(1,:), the_deg, eh(1,:))
set(gca, 'FontSize', 12, 'FontWeight', 'b')
xlabel('Incidence Angle \theta (deg)','FontSize',14, 'FontWeight','b')
ylabel('e', 'FontSize',14, 'FontWeight','b')
%title('Volume Scattering - Lossy Sand Particles', 'FontSize',16, 'FontWeight','b')
ch = get(gca, 'Children');
set(ch(2),'LineWidth', 1.5, 'LineStyle','-','Color','r')
set(ch(1),'LineWidth', 1.5, 'LineStyle','-','Color','b')
%axis([0 50 100 300])
legend ('v-pol', 'h-pol', 3)
grid



