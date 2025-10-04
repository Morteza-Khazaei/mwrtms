%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 5.1: sigma_note Polarization Synthesis of PolSAR Images

%Description: Given a polarimetric image, the code synthesizes the radar
%image for any specified transmit and receive polarization combination.

%Input Variables:
    % svv, shh, shv: 2-D complex arrays storing the vv-pol, hh-pol, and
    % hv-pol responses of the imaged scene. 
    %nx and ny: number of pixels of the image along x and y dimensions
    %psi_r: Rotation polarization angle of receive antenna of synthesized
           %image (degrees)
    %Chi_r: Ellipticity polarization angle of receive antenna of synthesized
            %image (degrees)
    %psi_t: Rotation polarization angle of transmit antenna of synthesized
            %image (degrees)
    %Chi_t: Ellipticity polarization angle of transmit antenna of synthesized
            %image (degrees)
            
%Output Product: 
    %Synthesized image.

%Book Reference: Section 5-10.3

%Matlab Code:

function [Im_out] = PolSynthesis_Images(svv, shh, shv, nx, ny, psi_r, chi_r, psi_t, chi_t)

Im_out = zeros(nx, ny);

% initialize variables
I_r = zeros(4,1);
I_t = zeros(4,1);
Q = zeros(4,4);

% transform angles to radian
d2r = pi/180;
two_psi_r = 2 *psi_r * d2r;
two_psi_t = 2* psi_t * d2r;
two_chi_r = 2* chi_r * d2r;
two_chi_t = 2* chi_t * d2r;

% calculate the Stokes vectors for the transmit and receive polarizations
I_r(1) = 0.5 *( 1 + cos(two_psi_r).* cos(two_chi_r));
I_r(2) = 0.5 *( 1 - cos(two_psi_r).* cos(two_chi_r));
I_r(3) = sin(two_psi_r) * cos(two_chi_r);
I_r(4) = sin(two_chi_r);

I_t(1) = 0.5 *( 1 + cos(two_psi_t).* cos(two_chi_t));
I_t(2) = 0.5 *( 1 - cos(two_psi_t).* cos(two_chi_t));
I_t(3) = sin(two_psi_t) * cos(two_chi_t);
I_t(4) = sin(two_chi_t);

% build transformation matrix Q
Q(1,1) = 1.0; Q(2,2) = 1.0;  Q(3,3) = 0.5; Q(4,4) = -0.5;

I_rQ = 4*pi * I_r' * Q; 

for ii = 1:nx
    for jj = 1:ny
        Pixel(2,2) = shh(ii,jj);
        Pixel(1,1) = svv(ii,jj);
        Pixel(1,2) = shv(ii,jj);
        Pixel(2,1) = Pixel(1,2);
        
        [sigma_note] = Sigma_note_Calc(I_rQ, I_t, Pixel);

        Im_out(ii,jj) = sigma_note;
    end
end

end

function [sigma_note] = Sigma_note_Calc(I_rQ, I_t, S)

%- S is a complex 2x2 matrix

%-- generate the modified Mueller matrix

M = zeros(4,4);
M(1,1) = abs(S(1,1)).^2;
M(1,2) = abs(S(1,2)).^2;
M(2,1) = abs(S(2,1)).^2;
M(2,2) = abs(S(2,2)).^2;

M(1,3) = real(conj(S(1,2))* S(1,1));
M(1,4) = -imag(conj(S(1,2))* S(1,1));

M(2,3) = real(conj(S(2,2))* S(2,1));
M(2,4) = -imag(conj(S(2,2))* S(2,1));

M(3,1) = 2* real(conj(S(2,1))* S(1,1));
M(4,1) = 2* imag(conj(S(2,1))* S(1,1));

M(3,2) = 2* real(conj(S(2,2))* S(1,2));
M(4,2) = 2* imag(conj(S(2,2))* S(1,2));

M(3,3) = real(S(1,1)*conj(S(2,2)) + S(1,2)*conj(S(2,1)));
M(4,4) = real(S(1,1)*conj(S(2,2)) - S(1,2)*conj(S(2,1)));

M(3,4) = -imag(S(1,1)*conj(S(2,2)) - S(1,2)*conj(S(2,1)));
M(4,3) = imag(S(1,1)*conj(S(2,2)) + S(1,2)*conj(S(2,1)));

%-- calculate sigma_note

sigma_note = I_rQ * M * I_t;

end

