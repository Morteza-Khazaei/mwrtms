function [vvdb,hvdb,vhdb,hhdb,ot,albedo,epsr_eff,...
    vv_vol,hv_vol,vh_vol,hh_vol,vv_surf,hv_surf,vh_surf,hh_surf] = ...
    DMRT_QMS_active(fGHz,deg0inc,depth,rho,dia,tau,...
    Tsnow,epsr_ground,rms,ratio,surf_model)
% active remote sensing of layered snowpack
% solve Active DMRT equation using quadrature eigenanalysis approach,
% assuming QCA Mie Sticky spheres (QMS)scattering model
% surface roughness is only accounted for surface backscattering.
% input arguments:
%   fGHz - frequency in GHz, scalar
%   deg0inc - incidence angle in degree, scalar
%   depth - layer thickness in centimeter, 1D array, from top to bottom
%   rho - snow density in gm/cc, 1D array
%   dia - ice grain diameter in centimeter, 1D array
%   tau - stickiness parameter in QCA, 1D array
%   Tsnow - snow temperature in (K), 1D array, (used to calc permittivity)
%   epsr_ground - ground relative permittivity
%   rms - rough ground rms height in cm, 0 assumes flat bottom boundary
%   ratio - correlation length divided by rms height
%   surf_model - surface scattering model, support physical model of 
%       'NMM3D' (through a look up table) and 1st order 'SPM3D'. 
%       Also support empirical 'OH' model.
%       * The range spanned by NMM3D LUT is 
%       rms height / wavelength in top medium:  (0.02101, 0.2101)
%       real part of permittivity of bottom media divided by top media: (3.0, 30.0)
%       ratio (correlation length / rms height): (7,15)
%       * OH model is validated for ks: 0.1~6, kl: 2.5~20, mv: 0.09~0.31
%       and incidence angle: 20~70 degrees
% 
% output parameters:
%   vvdb,hvdb,vhdb,hhdb - total backscatter in dB scale
%   
%   ot,albedo,epsr_eff: optical thickness, scattering albedo and effective
%       permittivity of each layer
%   vv_vol,hv_vol,vh_vol,hh_vol - volume scattering contribution to 
%       backscattering in linear scale
%   vv_surf,hv_surf,vh_surf,hh_surf - surface scattering contribution to 
%       backscattering in linear scale
% 
% Copyright: Electrical Engineering Department, University of Washington
% Revison Date:   July 28, 2014, Sep. 10, 2014
% 
% References: Tsang et al. TGRS 2007, Chang et al. JSTARS 2014,  
%   Liang et al., IGARSS 2008.
% 

Nquad = 32; % number of quadrature angles
ndeg = Nquad/2;

%% Parameters Setup
eps0 = 8.854e-12;
mu0 = pi*4e-7;
c = 1/sqrt(eps0*mu0);
wave = c/fGHz*1e-7; % cm
k = 2*pi/wave; % free space wave number in 1/cm

rho_ice = 0.917;

Nmax = fix(k*max(dia)) + 1;
Mmax = 4 * Nmax;

%% Calculate phase matrix from QCA
nl = length(depth);
np = 65;
P11_a = zeros(nl,np);
P22_a = zeros(nl,np);
P33_a = zeros(nl,np);
P34_a = zeros(nl,np);
P12_a = zeros(nl,np);
P21_a = zeros(nl,np);

perm_eff_a = zeros(nl,1);
ks_a = zeros(nl,1);
ka_a = zeros(nl,1);

for il = 1:nl
    epsr_ice = diel_ice(fGHz,Tsnow(il));
    fv = rho(il)/rho_ice;
    fv = max(0,min(1,fv));
    [keff,ka,ks,bmua,P11,P22,P33,P44,P34,P43] = PhaseMatrix12f_QMS(k,epsr_ice,dia(il),fv,tau(il),np);
    perm_eff_a(il)=(real(keff)/k)^2;  
    ka_a(il) = ka;
    ks_a(il) = ks;
    P11_a(il,:) = P11/ks;
    P22_a(il,:) = P22/ks;
    P33_a(il,:) = P33/ks;
    P34_a(il,:) = P34/ks;
end

ke_a = ks_a + ka_a;
albedo=ks_a./ke_a;
ot = ke_a.*depth(:);
epsr_eff = perm_eff_a;

%% surface backscattering
eps1 = perm_eff_a(end); % the bottom snow layer above ground
k1 = sqrt(eps1)*k;
k2=sqrt(epsr_ground)*k;
taiinc = asin(sind(deg0inc)./sqrt(perm_eff_a));  % incidence angle in snow
if rms > 0
    if strcmp(surf_model,'SPM3D')
        % surface scattering with 1st order SPM
        [svv,shh] = SPM3D(k1, taiinc(end), k2, rms/100, ratio);  % apply 1st order SPM with exponential correlation function
        sx = 0;
    elseif strcmp(surf_model,'NMM3D')
        % surface scattering from NMM3D LUT at 40 degree incidence angle
        [svv_dB, shh_dB, sx_dB] = NMM3D_LUT_NRCS_40degree_interp(rms/wave*sqrt(eps1),real(epsr_ground)/eps1,ratio);
        svv = 10^(svv_dB/10);
        shh = 10^(shh_dB/10);
        sx = 10^(sx_dB/10);
    elseif strcmp(surf_model,'OH')
        [svv,shh,sx] = rough_bks_OH(eps1,epsr_ground,taiinc(end),k1*rms);
    else
        error('Unsupported surface scattering model');
    end
else % flat surface
    svv = 0;
    shh = 0;
    sx = 0;
end

% fprintf('bare surface scattering: v:%.4f h:%.4f x:%.4f\n',svv,shh,sx);

% sigma_s = [svv sx; sx shh];
att = exp(-2*sum(ot./cos(taiinc)));
% sigma_s_a = sigma_s*att;
vv_surf = svv*att;
hh_surf = shh*att;
x_surf = sx*att;
hv_surf = x_surf;
vh_surf = x_surf;

%% calculate bistaticscattering pattern and backscattering
[sigma,gamma,deg0int] ...
    = AMLDRS(deg0inc, depth(:)', epsr_ground, [], ndeg, Mmax,...
    perm_eff_a,ka_a,ks_a,bmua,P11_a,P22_a,P33_a,P34_a,P12_a,P21_a); 

% volume scattering
sigma_ = 10.^(sigma/10);
vv_vol = sigma_(1);
hv_vol = sigma_(2);
vh_vol = sigma_(3);
hh_vol = sigma_(4);

% add surface scattering
vv = vv_vol + vv_surf;
hv = hv_vol + hv_surf;
hh = hh_vol + hh_surf;
vh = vh_vol + vh_surf;

vvdb=10*log10(vv);
hvdb=10*log10(hv);
hhdb=10*log10(hh);
vhdb=10*log10(vh);

% fprintf('vv = %.4f; (vol: %.4f, surf: %.4f)\t\t vvdb = %1.4f\n',vv,vv_vol,vv_surf,vvdb);
% fprintf('hv = %.4f; (vol: %.4f, surf: %.4f)\t\t hvdb = %1.4f\n',hv,hv_vol,x_surf,hvdb);
% fprintf('hh = %.4f; (vol: %.4f, surf: %.4f)\t\t hhdb = %1.4f\n',hh,hh_vol,hh_surf,hhdb);
% fprintf('vh = %.4f; (vol: %.4f, surf: %.4f)\t\t vhdb = %1.4f\n',vh,vh_vol,x_surf,vhdb);

% % bistatic scattering pattern
% gm_vv_db = gamma(1,:);
% gm_hv_db = gamma(2,:);
% gm_vh_db = gamma(3,:);
% gm_hh_db = gamma(4,:);
% 
% figure; subplot(2,1,1); plot(deg0int,gm_vv_db,'b','linewidth',2);
% hold on; plot(deg0int,gm_hh_db,'r--','linewidth',2);
% xlim([-70,70]); legend('VV','HH');
% subplot(2,1,2); plot(deg0int,gm_vh_db,'k--','linewidth',2);
% hold on; plot(deg0int,gm_hv_db,'g:','linewidth',2);
% xlim([-70,70]); legend('VH','HV');
% xlabel('scattering angle (degree)'); ylabel('Bistaticscattering coefficient (dB)'); 
end