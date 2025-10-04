% Demonstration of invoking DMRT_QMS_passive
addpath '..\common'
addpath '..'

% ==========================================================
% radiometer parameters
fGHz = 18.7; %36.5; % frequncy in GHz

% ==========================================================
% snowpack specification
% dia = 0.1;          % grain diameter in cm
% rho = 0.276;        % snow density in gm/cc
% tau = 0.1;          % QCA stikiness parameter
% depth = 30;         % snow layer thickness in cm
% Tsnow = 260;        % snow temperature in K

% OR load snowpack description file
[depth, rho, Tsnow, dia, tau] = load_snowpack('snowpack.txt');

% ==========================================================
% bottom boundary specification, temperature and permittivity
Tg = 270;           % ground temperature
% soil permittivity
mv = 0.02;
clayfrac = 0.3;
epsr_ground = soil_perm_MBSDM_Mironov(mv,clayfrac,fGHz);

% ==========================================================
% bottom boundary roughness specification
% * option 1: Q/H model
rough.model = 'QH';
rough.Q = 0;          % polarization mixing factor, unitless
rough.H = 0;          % roughness height factor, unitless
                      % Q = H = 0, means flat bottom surface        

% * option 2: Wegmuller and Matzler 1999 model
% rough.model = 'WM';
% rough.s = .000;       % rms height in meter

% ==========================================================
% invoke DMRT_QMS_passive to calculate brightness temperature
[TBv,TBh,deg0,ot,albedo,epsr_snow] = ....
    DMRT_QMS_passive(fGHz,dia,rho,tau,depth,Tsnow,Tg,epsr_ground,rough);

% ==========================================================
% observation
% Tb angular pattern
nangle = 91;
ob_anglei = linspace(0,90,nangle); % set up observation angles
ob_anglei = ob_anglei(ob_anglei <= deg0(end) & ob_anglei >= deg0(1));
Tb_v = spline(deg0,TBv,ob_anglei);
Tb_h = spline(deg0,TBh,ob_anglei);

figure; plot(ob_anglei,Tb_v,'-b','linewidth',2);
hold on; plot(ob_anglei,Tb_h,'--r','linewidth',2);
xlabel('observation angle (degree)');
ylabel('Brightness Temperature (K)');
legend('TB_V','TB_H'); xlim([ob_anglei(1) ob_anglei(end)])

% Tb in specifiend direction 
ob_angle = 55;  % observation angle in degree
Tb_v0 = spline(deg0,TBv,ob_angle);
Tb_h0 = spline(deg0,TBh,ob_angle);
fprintf('TB at %.1f degree: v %.2f(K); h %.2f(K)\n',ob_angle,Tb_v0,Tb_h0);