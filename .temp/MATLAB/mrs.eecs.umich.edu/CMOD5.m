%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
% Code: 16.1: The CMOD5 Geophysical Model Function

% Description: Code computes C-band VV $\sigma^0$ for a specified incidence 
% angle, radar azimuth angle, wind direction, and wind speed according to 
% the CMOD5 geophysical model function.
% 
% Input Variables:
% theta: incidence angle [15 <= theta <= 60] (deg)
% phi_radar: radar illumination azimuth angle relative to north (deg)
% phi_wind: wind direction relative to north (meteorological convention, deg)
% U: wind speed [0<U<60] (m/s)
% 
% Output Products:
% sigma^0:  VV normalized radar cross-section (unitless)
% 
% Book Reference: Chapter 16
% 

function sigma0=CMOD5(inarg)
% 
% sigma0=CMOD5(inarg)
%
% evaluate CMOD5 at particular wind speed & direction, azimuth angle,
% and incidence angle at the parameters specified by the rows of inarg
% where a row of inarg is given by 
%   inarg(1,:) = [U,D,A,T]
% with 
%   U is the wind speed in m/s
%   D is the wind direction relative to north in deg
%   A is the antenna azimuth angle relative to north in deg
%   T is the incidence angle in deg
%
% sigma-0 output is VV pol in linear space (not dB)
%

% extract columns of input rows
U=inarg(:,1)';  % speed
D=inarg(:,2)';  % wind direction
A=inarg(:,3)';  % radar azimuth angle
T=inarg(:,4)';  % incidence angle

% original CMOD5 constants
ZPOW = 1.6;
C= [-0.6880, -0.793, 0.338, -0.173, 0.00,  0.004, 0.111, ...
     0.0162,  6.340, 2.570, -2.180, 0.40, -0.600, 0.045, ...
     0.0070,  0.330, 0.012, 22.000, 1.95,  3.000, 8.390, ...
    -3.4400,  1.360, 5.350,  1.990, 0.29,  3.800, 1.530];

% evaluate model using matlab vectorized code
CChi     = cos((D - A)*pi/180);
X        = (T - 40.0) / 25.0; 
X2       = X.^2;
X3       = X.^3;
A0       = C( 1)+C( 2)*X+C( 3)*X2+C( 4)*X3;
A1       = C( 5)+C( 6)*X;
A2       = C( 7)+C( 8)*X;
GAM      = C( 9)+C(10)*X+C(11)*X2;
S0       = C(12)+C(13)*X;
S        = A2.*U;
S2       = S;
S2(S0>S) = S0(S0>S);
A4       = 1./(1+exp(-S2));
A3       = A4;
A3(S<=S0)=A4(S<=S0).*(S(S<=S0)./S0(S<=S0)).^(S0(S<=S0).*(1-A4(S<=S0)));
U0       = C(21) + C(22)*X + C(23)*X2;
D1       = C(24) + C(25)*X + C(26)*X2;
D2       = C(27) + C(28)*X;
A        = C(19)-(C(19)-1)/C(20);
B        = 1/(C(20)*(C(19)-1)^(C(20)-1));
B0       = (A3.^GAM).*10.^(A0+A1.*U);
B3       = C(15)*U.*(0.5+X-tanh(4*(X+C(16)+C(17)*U)));
B1       = (C(14)*(1+X)-B3)./(exp(0.34*(U-C(18)))+1);
U1       = U./U0+1;
U2       = U1;
U2(U1<C(19)) = A+B*(U1(U1<C(19))-1).^C(20);
B2       = (D2.*U2-D1).*exp(-U2);
sigma0   = B0.*(1+B1.*CChi+B2.*(2.*CChi.^2-1)).^ZPOW;
end
