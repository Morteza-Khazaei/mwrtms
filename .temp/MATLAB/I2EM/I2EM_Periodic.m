%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 10.4: I2EM Backscattering from Periodic Sinusoidal Surface

%Description: Code computes sigma_0_vv, hh, and hv, for a sinusoidal
%periodic surface with a random smaller-scale component.

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
    %sigma_0_vv(theta_0, phi_0), sigma_0_hh, sigma_0_hv
    
    %To maintain reasonable execution time for this code, calculation of
    %sigma_0_hv is switched off here. However, it can be activated with
    %ease (check the function called Calc_scatt_coeff listed below).
    
%Book Reference: Section 10-4.2

%Matlab code:

function [sigma_0_vv sigma_0_hh sigma_0_vh] = I2EM_Periodic(theta_0,phi_0 ...
    , eps, f, s, l, sp, Gmm, A)

%-- Note: the angles can be vectors or scalars

ntheta = length(theta_0);
nphi = length(phi_0);

sigma_0_vv = zeros(nphi, ntheta); sigma_0_hh = sigma_0_vv; sigma_0_vh = sigma_0_vv;


%-- compute the backscattering coefficient of the surface due to small
%scale roughness (no periodicity at this point)

theta_flat = 0:1:89; %local incidence angle in degrees (flat surface)
nflat = length(theta_flat);
theta_flat(1) = 0.01;

sig_f_vv = zeros(nflat, 1); sig_f_hh = sig_f_vv;

for n = 1: nflat
  [sig_f_vv(n) sig_f_hh(n)] = I2EM_Bistat_model(f, s, l, theta_flat(n), ...
      theta_flat(n), 180, eps, sp, 1);
end

sig_f_vv = 10.^(sig_f_vv/10); % transform to linear scale
sig_f_hh = 10.^(sig_f_hh/10);

t = 'finished calculating scattering of flat rough surface'

%-- set number of segments used to discretize the surface
ny = 100; % number of segments or points used to discretize the surface
y = linspace(0, Gmm, ny);

Zx = 0; %partial derivative of periodic surface along x

sig_ss_vv = zeros(ny,1); 
sig_ss_hh = sig_ss_vv; sig_ss_vh = sig_ss_vv;


%--- calculate the scattered response from the periodic surface

for nt = 1 : ntheta
    for np = 1: nphi
 
        for ii = 1: ny
%            Zy = 2*pi/Gmm *A * sin(2*pi/Gmm * y(ii));
            Zy = 2*pi/Gmm *A * cos(2*pi/Gmm * y(ii));
            alpha = atan(Zy);
    
           [sig_ss_vv(ii), sig_ss_hh(ii), sig_ss_vh(ii)] = Calc_scatt_coeff(Zx,Zy, ...
        theta_0(nt), phi_0(np), sig_f_vv, sig_f_hh, theta_flat);
     
           sig_ss_vv(ii) = sig_ss_vv(ii) * sec(alpha); 
           sig_ss_hh(ii) = sig_ss_hh(ii) * sec(alpha); 
           sig_ss_vh(ii) = sig_ss_vh(ii) * sec(alpha); 
        end

        sigma_0_vv(np,nt) = trapz(y,sig_ss_vv) / Gmm; % integrate
        sigma_0_hh(np,nt) = trapz(y,sig_ss_hh) / Gmm;
        sigma_0_vh(np,nt) = trapz(y,sig_ss_vh) / Gmm;
    end
end

end

%-------------------------------------------------------------------------

function [sig_ss_vv, sig_ss_hh, sig_ss_vh] = Calc_scatt_coeff(Zx,Zy,th,ph,...
    sig_f_vv, sig_f_hh, theta_flat)

%-It calculates the scattering coefficients of the differential area at a
%-given point on the periodic surface due to rough surface scattering from
%-the small scale roughness. 


%-calculate the polarization vectors in both local and global coordinates
th_rad = th*pi/180; %present angle in radian
ph_rad = ph*pi/180;

[v,h,v_pr,h_pr,thpr] = vectors2(Zx,Zy,th_rad, ph_rad);

%-calculate local backscattering coefficients due to small scale roughness

  % --- this part is not used anymore. It still does work ----%
% % % [sig_s_vv, sig_s_hh, sig_s_vvhh, sig_s_hv] = Calc_localScat_PO(thpr, conj(eps), f, s, l);
% % % 
% % % sig_ss_vv = dot(v, v_pr).^4 .* sig_s_vv + dot(v, h_pr).^4 .*sig_s_hh + ...
% % %     2 * dot(v, h_pr).^2 .* dot(v, v_pr).^2 .* sig_s_vvhh;
% % % 
% % % sig_ss_hh = dot(h, v_pr).^4 .* sig_s_vv + dot(h, h_pr).^4 .*sig_s_hh + ...
% % %     2 * dot(h, h_pr).^2 .* dot(h, v_pr).^2 .* sig_s_vvhh;
% % % 
% % % sig_ss_vh = dot(v, v_pr).^2 .* dot(v_pr, h).^2 .* sig_s_vv + dot(v, h_pr).^2 ...
% % %     .* dot(h_pr, h).^2 .*sig_s_hh + 2 * dot(v, v_pr) .* dot(h, h_pr) ... 
% % %     .* dot(h, v_pr) .* dot(h_pr, v) .* sig_s_vvhh + ( dot(v, v_pr).* ...
% % %     dot(h_pr, h) + dot(v, h_pr) .* dot(v_pr, h)).^2 .* sig_s_hv;
  %-----------------------------------------------------------%

thpr_deg = thpr*180/pi; % need this angle in degrees

if thpr_deg > 89

    sig_s_vv = 0; %No return from this portion of the surface 
    sig_s_hh = 0;
else    
    if thpr_deg < theta_flat(1)
        thpr_deg = theta_flat(1);
    end
%--extract the correct data from calculated scattering of flat surface
    sig_s_vv = interp1(theta_flat, sig_f_vv, thpr_deg); 
    sig_s_hh = interp1(theta_flat, sig_f_hh, thpr_deg); 

end
  %-- Calculation of x-pol can be done as shown below. However, it will
  %slow down the entire code considerably. It have been commented out here.
  %However, the user may choose to use it by removing the comments.
%--calculate x-pol using IEMX_model. 
%  [sig_s_hv] = IEMX_model(f, s, l, thpr_deg, eps, sp,1, 0);

sig_s_hv = 10e-20; % We are discontinuing calc of x-pol 

%-- The following expressions keep only the strong terms.
% sig_ss_vv = dot(v, v_pr).^4 .* sig_s_vv ;
% sig_ss_hh = dot(h, h_pr).^4 .*sig_s_hh ;

sig_ss_vv = dot(v, v_pr).^4 .* sig_s_vv + dot(v, h_pr).^4 .*sig_s_hh  ...
   + 2 * dot(v, h_pr).^2 .* dot(v, v_pr).^2 .* sqrt(sig_s_vv * sig_s_hh);
sig_ss_hh = dot(h, h_pr).^4 .*sig_s_hh + dot(h, v_pr).^4 .* sig_s_vv ...
   + 2 * dot(h, h_pr).^2 .* dot(h, v_pr).^2 .* sqrt(sig_s_vv * sig_s_hh);

sig_ss_vh = (dot(v, v_pr).* dot(h_pr, h) + ...
            dot(v, h_pr) .* dot(v_pr, h)).^2 .* sig_s_hv;

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

% This function is no longer in use now. 
function [sig_s_vv, sig_s_hh, sig_s_vvhh, sig_s_hv] = Calc_localScat_PO(th, eps, f, s, l)

%-calculates the noncoherent backscattering coefficients of the small scale roughness
%surface based on Physical Optics model

[rvv, rhh]=refl_coef(th, 1, eps); %calculates reflection coefficients at selected angle

wn = 2*pi * f / 0.3; % wavenumber

%-calculate coherent component at th = 0 deg

d0 = 0.01; % width of the illuminated area

th_err = 0.005 * pi/180;

if th <=th_err
    sig_s_vv_c = (4*wn^2*d0^2 * abs(rvv)^2 /pi) .*exp(-(2*wn*s*cos(th))^2);
    sig_s_hh_c = (4*wn^2*d0^2 * abs(rhh)^2 /pi) .*exp(-(2*wn*s*cos(th))^2);
else
    sig_s_vv_c = 0;
    sig_s_hh_c = 0;
end

%--calculate incoherent component
xx = (2 * wn * s* cos(th) )^2;
sum = 0;

for ii = 1:31
    sum = sum + (xx)^ii ./ factorial(ii);
end

cons = 2*wn * sin(th);
intg = 0;

% delx = 0.0001;
% for ii = 1: 100000
%     x = (ii-1) * delx;
%     rho = exp(-x^2/(l1^4 + x^2*l2^2)^0.5);
%     ss = rho*besselj(0,x*cons)*x*delx;
%     intg = intg + ss;
%  %   x, rho, ss
% end
% intg = intg * sum;



% load smallscale_corrfunc.mat;
% b = max(xmod); 
% a = 0;
% xi = linspace(a,b,10000);
% rho = interp1(xmod, tt_r, xi);

b = 1; 
a = 0;
xi = linspace(a,b,10000);
% rho = exp(-xi.^2./(l1.^4 + xi.^2 .*l2.^2).^0.5); %gaussian/exponential
% rho = exp(-xi.^2 / l1.^2); %gaussian
% rho = exp(-xi / l1); %Exponential

 rho = (1 + xi.^2 ./l.^2) .^(-1.5); % ulaby's


[JJ, ierr] = besselj(0, xi .*cons); 
yi = rho.* JJ .*xi;

if ierr ~= 0
    ierr
end

%I1 = quad(@(x) interp1(xi, yi, x, 'cubic',0), a, b, 1.0e-12);

I1 = trapz(xi, yi);

intg = I1 .*sum;


%--calculate the backscattering coefficients

sig_s_vv_n  = 2 *( wn*abs(rvv)*cos(th))^2 * exp(-(2*wn*s*cos(th))^2) * intg;
sig_s_hh_n  = 2 *( wn*abs(rhh)*cos(th))^2 * exp(-(2*wn*s*cos(th))^2) * intg;


sig_s_vv = sig_s_vv_c + sig_s_vv_n;
sig_s_hh = sig_s_hh_c + sig_s_hh_n;


sig_s_hv = 0;

sig_s_vvhh = sig_s_vv * abs(real(rvv * conj(rhh)))/abs(rvv)^2;

end
function [rho_v, rho_h]=refl_coef(the1, eps1, eps2)

% calculates the v and h-polarized reflection coefficients of a plane
% dielectric surface

n1 = sqrt(eps1);
n2 = sqrt(eps2);
costh2 = sqrt(1 - (n1.*sin(the1)./n2).^2);


rho_v = -(n2.*cos(the1) - n1.*costh2)./(n2.*cos(the1)+n1.*costh2);
rho_h = (n1.*cos(the1) - n2.*costh2)./(n1.*cos(the1)+n2.*costh2);
end
