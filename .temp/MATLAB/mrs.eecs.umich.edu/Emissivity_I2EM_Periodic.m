%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 12.4: Emissivity of a periodic Sinusoidal Surface

%Description: Code computes emissivity for a random surface superimposed on
%a periodic  sinusoidal surface.

%Input Variables:
    %eps: dielectric constant of the scattering medium (complex)
    %sp: correlation function of the random component:
            % 1- Exponential, 2- Gaussian
    %s: rms height of the random component (m)
    %l: correlation length of the random component (m)
    %Gmm: Spatial period of the periodic surface (m)
    %A: amplitude of periodic sinusoid (m)
    %f: frequency  in GHz
    %theta_0: Incidence angle (deg)
    %phi_0: Azimuth angle(Fig. 10-24) (deg)

%Output Products:
    %e_v(theta_0, phi_0): v-polarized emissivity
    %e_h(theta_0, phi_0): h-polarized emissivity    
   
%Book Reference: Section 12-4

%Matlab code:

function [ev eh] = Emissivity_I2EM_Periodic(theta_0,phi_0 ...
    , eps, f, s, l, sp, Gmm, A)

ny = 100; % number of segments or points used to discretize the surface
y = linspace(0, Gmm, ny);

Zx = 0; %partial derivative of periodic surface along x

ev_ss = zeros(ny,1); 
eh_ss = ev_ss;

for ii = 1: ny
    Zy = 2*pi/Gmm *A * sin(2*pi/Gmm * y(ii));
    alpha = atan(Zy);
    
    [ev_ss(ii), eh_ss(ii)] = Calc_scatt_coeff(Zx,Zy, ...
        theta_0, phi_0, eps, f, s, l, sp);
     
     ev_ss(ii) = ev_ss(ii) * sec(alpha); 
     eh_ss(ii) = eh_ss(ii) * sec(alpha); 
end

ev = trapz(y,ev_ss) / Gmm; % integrate
eh = trapz(y,eh_ss) / Gmm; % integrate


end
%-------------------------------------------------------------------------

function [ev_ss eh_ss] = Calc_scatt_coeff(Zx,Zy,th,ph,eps, f, s, l, sp)

%-It calculates the emissivity of the differential area at a
%-given point on the periodic surface due to rough surface scattering from
%-the small scale roughness. 


%-calculate the polarization vectors in both local and global coordinates
th_rad = th*pi/180; %present angle in radian
ph_rad = ph*pi/180;

[v,h,v_pr,h_pr,thpr] = vectors2(Zx,Zy,th_rad, ph_rad);

%-calculate local emissivity  due to small scale roughness

thpr_deg = thpr*180/pi; % need this angle in degrees

%--calculate local emissivity due to small-scale roughness using I2EM
[ev_s eh_s] = Calc_emissivity_I2EMmodel(f, s, l, thpr_deg, eps, sp);


ev_ss = dot(v, h_pr).^2 .* eh_s + dot(v, v_pr).^2 .* ev_s ;
eh_ss = dot(h, h_pr).^2 .* eh_s + dot(h, v_pr).^2 .* ev_s ;


end

function [v,h,v_pr,h_pr,thpr] = vectors2(Zx,Zy,th, ph)
%
% Calculates the local coordinates and polarization vectors in global
% coordinate system. 

snt = sin(th);
cst = cos(th);
snp = sin(ph);
csp = cos(ph);

D0 = (1 + Zx.^2 + Zy.^2).^(-0.50);

%-- calculate polarization vectors, unit normal vector, and incident wave
%directions in global coordinates

n(1) = -Zx; % unit normal vector at point on periodic surface
n(2) = -Zy;
n(3) = 1;
n = n * D0;

n = n / norm(n);

h(1) = -snp;
h(2) = csp;
h(3) = 0;

v(1) = -cst.*csp ;
v(2) = -cst.*snp;
v(3) = -snt;

ni(1) = snt .* csp; %incident wave direction
ni(2) = snt .* snp;
ni(3) = -cst;
ni = ni / norm(ni);


costhpr = - dot(ni, n);
thpr = acos(costhpr);
sinthpr = sin(thpr);

%--calculate local coordinate axes
z_pr = n;

tmp = cross(n, ni);
y_pr = tmp /norm(tmp);

x_pr = cross(y_pr, z_pr);

h_pr = y_pr;

tmp = cross(h_pr, ni);
v_pr = tmp / norm(tmp);


end
