%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 3.2: Array factor for uniformly spaced array with linear phase
%distribution
%Description: Code computes the array factor as defined by eq. 3.97 for an
%N-element uniformly spaced array.
    
%Input Variables:
    %N: number of array elements
    %a: amplitude of feeding coefficients
    %d: spacing between adjacent elements (m)
    %f: frequency in GHz
    %delta: Incremental phase delay (degrees)
    
%Output Products:
    % normalized array factor Fa versus angle theta
    
%Book Reference: Section 3-9

%Example call: [Fa thetad] = ArrayFactor_LinearPhaseDist(N, a,d, f, delta)

%MATLAB Code

function [Fa thetad] = ArrayFactor_LinearPhaseDist(N, a,d, f, delta)

  delta = delta .*pi/180; % transform the phase delay to radian
  
  thetad = -180:1:180; % angle in degrees
  theta = thetad .*pi/180; % transform to radians
  
  k =  20*pi/3*f; % wave number 
    
  na = length(thetad); % length of angle vector
  Fa = zeros(na,1); %creat an array of same length
       
  for ii = 1: na
      fact = k * d* cos(theta(ii)) - delta;
      
      for jj = 1: N          
         Fa(ii) = Fa(ii) + a(jj) .* exp(1i*(jj-1).*fact);
      end
  end
  Fa = (abs(Fa)).^2;
  
  end