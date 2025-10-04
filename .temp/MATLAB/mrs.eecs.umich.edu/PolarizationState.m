%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 2.1: Polarization State 
%Description: Code computes the angles of polarization ellipse of EM plane
%wave traveling in +z direction.
%Input Variables:
    %a_x: amplitude of x-component of E-field (arbitrary units)
    %a_y: amplitude of y-component of E-field (arbitrary units)
    %delta: phase of y-component relative to phase of x-component (phi_y -phi_x) (in degrees)

%Output Products:
    %psi: rotation angle of polarization ellipse (degrees)
    %chi: ellipticity angle of polarization ellipse (degrees)
%Book Reference: Section 2-3
%MATLAB Code: PolarizationState.m

%Example call: [psi chi] = PolarizationState(a_x,a_y,delta)

%MATAB CODE

function [psi chi] = PolarizationState(a_x, a_y, delta )

delta = delta * pi/180; % transform the phase angle to radian

alpha = atan2(a_y, a_x); % calculate the auxiliary angle (0<alpha<pi/)

cosdelta = cos(delta);
    
twopsi = abs(atan( tan(2* alpha) .* cos(delta))); 

if alpha > pi/4
    psi = 0.5 * (pi - twopsi);
else
    psi = 0.5 * twopsi;
end

chi = 0.5 * asin( sin(2*alpha) .* sin(delta)); 


if cosdelta > 0
    psi = abs(psi);
end
if cosdelta < 0
    psi = - abs(psi);
end

psi = psi * 180/pi; % transform to degrees
chi = chi * 180/pi; % transform to degrees

end

