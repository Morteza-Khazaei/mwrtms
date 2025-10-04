%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

f = 1:.01:1000; % freq in GHz
Rr = 2.54; % rainfall rate in mm/hr

nf = length(f);

extinc = zeros(nf,1);

for ii = 1: nf
   extinc(ii) = RainExtincCoef(Rr,f(ii));
end

figure(1)
loglog(f, extinc)
axis([ 1 1000 0.001 100])
grid

%---------------------------------------

f = 35; 
Rr = 0.5:0.5:100;

nf = length(Rr);

extinc = zeros(nf,1);

for ii = 1: nf
   extinc(ii) = RainExtincCoef(Rr(ii),f);
end

figure(2)
plot(Rr, extinc)
grid
