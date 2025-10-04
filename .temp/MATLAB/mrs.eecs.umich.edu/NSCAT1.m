%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
% Code: 16.2: The NSCAT1 Geophysical Model Function

% Description: Code computes Ku-band VV or HH $\sigma^0$ for a specified 
% incidence angle, radar azimuth angle, wind direction, wind speed, and 
% polarization according to the NSCAT1 geophysical model function.
% 
% Input Variables:
% theta: incidence angle [15 <= theta <= 60] (deg)
% phi_radar: radar illumination azimuth angle relative to north (deg)
% phi_wind: wind direction relative to north (oceanographic convention, deg)
% U: wind speed [1<= U <= 30](m/s)
% p: polarization flag (HH: p=1;  VV: p=2)
% 
% Output Products:
% \sigma^0:  normalized radar cross-section (unitless)
% 
% Book Reference: Chapter 16 


function sigma0=NSCAT1(inarg)
% 
% sigma0=NSCAT1(inarg)
%
% Evaluate NSCAT1 at particular wind speed & direction, azimuth angle,
% incidence angle, and polarization  at the parameters specified by 
% the rows of inarg where a row of inarg is given by 
%   inarg(1,:) = [U,D,A,T,p]
% with 
%   U is the wind speed in m/s (1-30)
%   D is the wind direction relative to north in deg
%   A is the antenna azimuth angle relative to north in deg
%   T is the incidence angle in deg (16-66)
%   p is the polarization (1=vv, 2=hh)
%
% returned sigma-0 output is in linear space (not dB)

% this routine reads the file nscat1.txt which contains a table
% of Ku-band VV and HH sigma0 values

% to save time for multiple calls, check to see if tabular model 
% function has been previously loaded into global memory
if exist('nscat1_gmf_model_tab','var')==0
  %  load tabular model function into global memory
  file='nscat1.txt';
  tmp=load(file);
  global NSCAT1_GMF
  NSCAT1_GMF=reshape(tmp,[26,50,37,2]);
else
  % enable access to previously loaded table
  global NSCAT1_GMF
end

% extract parameters from input array
speed=inarg(:,1);
chi=mod(360+inarg(:,2)-inarg(:,3),360);
inc=inarg(:,4);
ipol=inarg(:,5);
sigma0=zeros(size(speed));

% polarization
kpol=ipol*0+1;
kpol(ipol==2)=2;

% adjust the relative azimuth into the range (0..180) degrees
chi(chi<0)=chi(chi<0)+360;
chi(chi>180)=360-chi(chi>180);
ichi=floor(chi/5);
q=chi/5-ichi;
q1=1-q;

% interpolation coefficients for interpolation in incidence angle
htht = inc/2;
kt1  = floor(htht);
kt1(kt1< 8) = 8;
kt1(kt1>32) = 32;
kt2 = kt1 + 1;
at2 = htht - kt1;
at1 = 1 - at2;
kt1 = kt1-7;
kt2 = kt2-7;

% interpolation coefficients for interpolation of tabular speed
speed(speed<1)=1;
speed(speed>30)=30;
ju =floor(speed);
p = speed - ju;
p1= 1 - p;

% multi-linear interpolation of the tabular GMF
mfun=zeros([2,2]);
for k=1:length(speed)
  for i = 1:2
    kchi = i + ichi(k) - 1;
    kchi(kchi>36)=72-kchi(kchi>36);
    kchi(kchi<1)=1;
    for j = 1:2
      ku = j + ju(k) - 1;
      if ku==0
	mfun(j,i)=0;
      else
	mfun(j,i) = at1(k)*NSCAT1_GMF(kt1(k),ku,kchi,kpol(k)) + ...
	            at2(k)*NSCAT1_GMF(kt2(k),ku,kchi,kpol(k));
      end
    end
  end

  sigma0(k) = p1(k)*q1(k)*mfun(1,1) + p(k)*q1(k)*mfun(2,1) + ...
               q(k)*p1(k)*mfun(1,2) + p(k)* q(k)*mfun(2,2);
end
