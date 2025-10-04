%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

% To use this code without significant modifications, the user is required
% to provide a file storing the polarimetic response of a scene (i.e. complex
% data). A total of nx.ny rows of complex data should be stored in the
% file, where nx and ny correspond to the number of pixels in x and y
% dimensions, respectively. In each row the data stored is as follows: 
    % Re(Shh) Im(Shh) Re(Shv) Im(Shv) Re(Svv) Im(Svv)
%
% Replace the file name aa.dat with the name of the new data set and modify
% the dimensions of the image accordingly. 
%
% You can then modify the polarimetric angles (see 2nd evaluation cell below) 
% that characterize the transmit and receive polarizations of the
% synthesized image. 


%-- load the original file holding the polarimetric image data

load aa.dat;
mydata = aa;

nx = 512; ny = 512; % dimensions of image

svv = zeros(nx,ny); % initialize arrays
shh = svv; shv = svv;

for ii = 1:nx % Cast the data into nx x ny arrays and store the scattering 
              %matrix elements into separate arrays.
    for jj = 1:ny
        na = jj + ny *(ii-1);
        svv(ii,jj) = mydata(na,5) + 1i* mydata(na,6);
        shv(ii,jj) = mydata(na,3) + 1i* mydata(na,4);
        shh(ii,jj) = mydata(na,1) + 1i* mydata(na,2);
    end
end

orig_image = 10*log10(abs(svv).^2); % in dB scale

max1 = max(max(orig_image)); % maximum return from image
orig_image = orig_image - max1; %normalize by maximum value

figure(1)
imagesc(orig_image, [-40 0])
title('Original Image: VV-polarization')
colorbar


%% 2nd evaluation cell

%---------Begin constructing the synthesized image here ------------------
%-- Transmit and receive polarization angles of the synthesized image

%- transmit:
psi_t = 90; 
chi_t = 45;

%- receive:
psi_r = 0;
chi_r = 45;


[new_image] = PolSynthesis_Images(svv, shh, shv, nx, ny, psi_r, chi_r, psi_t, chi_t);

new_image = 10*log10(new_image); % cast in dB scale

max1 = max(max(new_image)); % maximum return from image
new_image = new_image - max1; %normalize by maximum value

figure(2)
imagesc(new_image, [-40 0])
title('Synthesized Image')
colorbar
