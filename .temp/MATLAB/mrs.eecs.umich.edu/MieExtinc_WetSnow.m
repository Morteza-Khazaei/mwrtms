%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 11.4: Mie Extinction of Wet Snow

%Description: Code computes the Mie absorption, scattering, and extinction
%coefficients of wet snow for specified snow density and ice particle
%radius. The wet snow temperature is assumed to be 0 degrees C.

%Input Variables:
    % rho_s: Snow density (g/cm3)
    % ri: Ice-particle radius (m)
    % f: frequency (GHz)
    % mv: snow wetness (%)
    
%Output Products:
    % kappa_a: absorption coefficient (Np/m)
    % kappa_s: scattering coefficient (Np/m)
    % kappa_e: extinction coefficient (Np/m)
    % a: single-scattering albedo
    
%Book Reference: Section 11-15.1 and eq 11.112 and 11.113 with Q computed
%according to the Mie model of section 8-5.

%Matlab code

function [kappa_a kappa_s kappa_e a] = MieExtinc_WetSnow(rho_s,ri,f,mv)

mv= mv ./100; % transform from percentage to absolute.
t = 0; %set temperature to 0 degree C.

rho_i = 0.9167; %density of ice (g /cm3)

%- calculate relative dielectric constant of pure ice
[epsr epsi] = RelDielConst_PureIce(t,f);
eps_sp = epsr -1i * epsi;

%- calculate the relative dielectric constant of water 
[epsr, epsi]= RelDielConst_PureWater(t,f);
eps_w = epsr -1i*epsi; %cast into complex 

%- Calculating the effective dielectric constant of the wet-air medium
%(i.e. the background medium of wet snow):

%according to the mixing formula of water/air mixture (eq: 11.114), and
%assuming that the depolarization factors are given as Aa = Ab = 0.06, and
%Ac = 0.88, then the mixing formula of the effective dielectric constant of
%the background medium, eps_b, reduces to third order polynomial: 
%   a *eps_b^3 + b*eps_b^2 + c*eps_b + d = 0:

Aa = 0.06; Ac = 0.88;

c3 = (1 - Ac - Aa + Ac * Aa); %(1-Aa)*(1-Ac);
c2 = eps_w * (Ac + Aa - 2* Ac*Aa); %eps_w*( (1-Ac)*Aa+(1-Aa)*Ac ); 
c1 = Aa * Ac * eps_w.^2;
a1 = eps_w *( 2*Ac + Aa); %(2*Ac*+ Aa)*eps_w;
a2 = (3 - 2*Ac - Aa); %(2*(1-Ac)+(1-Aa));
B = mv/3*(eps_w-1);

%-next we cast the polynomial coefficient into a vector (X) and use the matlab
%function roots to retrieve the 3 different roots of the equation

X(1) = c3; % note these coefficients are complex numbers.
X(2) = c2 - c3 - B * a2;
X(3) = c1 - c2 - B*a1;
X(4) = -c1;

rt_b = roots(X); % roots of the cubic polynomial (with a <>0)

%- which root to select?
%-now we need to discard the roots with negative real parts. Set the root to
%very large quantity. Then we select the roots with the smallest absolute
%value.

% the following is simple testing to verify that the roots are indeed
% correct values.
% x = -rt_b(1) + 1 + mv/3.*(eps_w-1).* ( 2./(1 + Aa*(eps_w./rt_b(1)-1)) + ...
%     1./(1 + Ac.*(eps_w./rt_b(1)-1)))
% 
% x = -rt_b(2) + 1 + mv/3.*(eps_w-1).* ( 2./(1 + Aa*(eps_w./rt_b(2)-1)) + ...
%     1./(1 + Ac.*(eps_w./rt_b(2)-1)))
% 
% x = -rt_b(3) + 1 + mv/3.*(eps_w-1).* ( 2./(1 + Aa*(eps_w./rt_b(3)-1)) + ...
%     1./(1 + Ac.*(eps_w./rt_b(3)-1)))
% 
% 
% 
% y = X(1)*rt_b(1)^3 + X(2)*rt_b(1)^2 + X(3)*rt_b(1) + X(4)
% 
% y = X(1)*rt_b(2)^3 + X(2)*rt_b(2)^2 + X(3)*rt_b(2) + X(4)
% 
% y = X(1)*rt_b(3)^3 + X(2)*rt_b(3)^2 + X(3)*rt_b(3) + X(4)


nr = length(rt_b);

for ii = 1:nr
    if real(rt_b(ii)) < 0,
        rt_b(ii) = 1.0e6;
    end
end
rt_min = min(abs(rt_b));
eps_b = rt_b(rt_min == abs(rt_b));  % dielectric constant of background


%  eps_b1 = mv/3 *(eps_w - 1) + 1; %crude approximation
% 
% eps_b2 = TVBmodel_HeterogeneousMix(eps_w, 1, 1, mv); %Using TVB model
% (approximate solution)
% 
% eps_b3 = eps_w*(1+ 2/3*mv*(eps_w -1)) ./(eps_w-mv/3*(eps_w-1)); %using
% the disk approximation of water droplet shape in Polder Van Santen
% formula 
%

%-- calculate the Mie efficiencies using Code 8.12
% [Es Ea Ee Eb t1 t2 t3 t4] = ...
%    Mie_Rayleigh_ScatteringOfSpheres(ri, f, eps_sp, eps_b); 

[Es Ea Ee Eb t1 t2 t3 t4] = ...
   Mie_Rayleigh_ScatteringOfSpheres(ri, f, eps_sp, real(eps_b)); 

area = pi * ri^2;

Qs = Es * area;
Qa = Ea * area;
Qe = Ee * area;

vi = rho_s ./rho_i;
Nv = vi ./(4/3 * pi .*ri^3); 

%coefficients due to the scatteres:
kappa_ai = Nv .* Qa;
kappa_s = Nv .* Qs; 
% kappa_e = Nv .* Qe;

%-absorption by the background
kappa_ab = 2* (2*pi*f/0.3) .*(1-vi) .* abs(imag(sqrt(eps_b)));


kappa_a = kappa_ai + kappa_ab; %total absorption

kappa_e = kappa_s + kappa_a; % total extinction

a = kappa_s /kappa_e; %albedo

end

    