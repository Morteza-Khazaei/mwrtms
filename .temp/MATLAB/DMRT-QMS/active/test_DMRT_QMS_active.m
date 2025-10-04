% Demonstration of invoking DMRT_QMS_active
addpath '..\common'
addpath '..'

% ==========================================================
% SAR or scatterometer parameters
deg0inc = 40;                       % incidence angle in degree
fGHz = 17.2; % 9.6; 17.2; % 13.4    % radar frequency in GHz

% ==========================================================
% snowpack description
% depth = 20;         % snow depth in centimeter
% rho = 0.345;        % snow density in gm/cc
% dia = 1.5 / 10;     % ice grain diameter in centimeter
% tau = 0.1;          % QCA stickiness parameter
% Tsnow = 260;        % specify snow temperature to calculate permittivity

% OR load snowpack description file
[depth, rho, Tsnow, dia, tau] = load_snowpack('snowpack.txt');

% ==========================================================
% bottom boundary specification, permittivity
% soil permittivity
mv = 0.15; % 0.02;
clayfrac = 0.3;
epsr_ground = soil_perm_MBSDM_Mironov(mv,clayfrac,fGHz);
% epsr_ground = 5 + 0.5i;

% ==========================================================
% bottom boundary roughness specification
rms = 0.25;                 % rough ground rms height, (cm)
                            % rms == 0 assumes flat bottom boundary
ratio = 7;                  % correlation length / rms height

% physical models to calculate surface backscattering
% * option 1: 'NMM3D', (through a look up table)
%   The range spanned by NMM3D LUT is 
%   rms height / wavelength in top medium:  (0.02101, 0.2101)
%   permittivity of bottom media divided by top media: (3.0, 30.0)
%   ratio (correlation length / rms height): (4,15) 
surf_model = 'NMM3D';       % pre-built NMM3D look up table

% * option 2: 1st order 'SPM3D'
% surf_model = 'SPM3D';     % first order SPM

% empirical models to predict surface backscattering
% * option 3: 'OH' model 
% surf_model = 'OH';        % Oh and Sarabandi 1992, ratio not used  

% ==========================================================
% invoke DMRT_QMS_active to calculate backscatter
[vvdb,hvdb,vhdb,hhdb,ot,albedo,epsr_eff,...
    vv_vol,hv_vol,vh_vol,hh_vol,vv_surf,hv_surf,vh_surf,hh_surf] = ...
    DMRT_QMS_active(fGHz,deg0inc,depth,rho,dia,tau,...
    Tsnow,epsr_ground,rms,ratio,surf_model);

% outputs
fprintf('vv = %.4f; (vol: %.4f, surf: %.4f)\n',vvdb,10*log10(vv_vol),10*log10(vv_surf));
fprintf('hh = %.4f; (vol: %.4f, surf: %.4f)\n',hhdb,10*log10(hh_vol),10*log10(hh_surf));
fprintf('hv = %.4f; (vol: %.4f, surf: %.4f)\n',hvdb,10*log10(hv_vol),10*log10(hv_surf));
fprintf('vh = %.4f; (vol: %.4f, surf: %.4f)\n',vhdb,10*log10(vh_vol),10*log10(vh_surf));