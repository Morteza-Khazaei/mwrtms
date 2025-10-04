%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 4.6: Relative Dielectric Constant of Wet Snow
%Description: Code computes the real and imaginary parts of the relative
%dielectric constant of Wet Snow
%Input Variables:
    %Ps: Dry Snow Density (g/cm^3)
    %mv: volumetric water content (0<mv<30 in percentage scale)
    %f: frequency in GHz
%Output Products:
    %epsr: real part of relative permitivity
    %epsi: imaginary part of relative permitivity
%Book Reference: Section 4-6.2
%MATLAB Code: RelDielConst_WetSnow.m

%Example call: [A B] = RelDielConst_WetSnow(Ps,mv, f)
%Computes the real and imaginary components of the permitivity of Wet Snow
    %based on the dry snow density (Ps), volumetric water content (mv), and frequency
    %vector (f) and assigns them to vectors A and B respectively.

%MATAB CODE

function [epsr epsi] = RelDielConst_WetSnow(Ps, mv, f )

A1 =  0.78 + 0.03 .*f - 0.58e-3 .* f.^2;
A2 =  0.97 - 0.39e-2 .*f + 0.39e-3 .* f.^2;
B1 =  0.31 - 0.05 .*f + 0.87e-3 .* f.^2;

A = A1 .*(1.0 + 1.83*Ps + 0.02*mv^1.015) + B1;
B = 0.073 .*A1;
C = 0.073 .*A2;
x = 1.31;
f0 = 9.07;

epsr = A + (B .* mv^x) ./ (1+(f/f0).^2);
epsi = (C .* (f/f0) .* mv^x) ./ (1+(f/f0).^2);
end
