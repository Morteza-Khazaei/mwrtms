%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

% note that in Matlab the array or vector index cannot start with 0, so the
% the vectors in this code are constructed with the first element indexed
% with 1. The equations reported in the text and used in this code will be
% modified accordingly. 

 f = 10; % frequency in GHz
 N = 15; %Number of array elements

 d = 0.015; % spacing between adjacent elements in m

 delta = 50;  % phase delay in degrees 
 
%--case 1: equal amplitude array 
a = ones(N,1); % amplitude set to 1 


[Fa1 theta1] = ArrayFactor_LinearPhaseDist(N,a,d,f, delta);
Fa1 = Fa1 ./max(Fa1); % normalize the array

%--case 2: binomial distribution

for ii = 1:(N+1)/2
    a(ii) = factorial(N-1) ./( factorial(ii-1) .* factorial(N-ii-2));
    a(N-ii) = a(ii);
end

[Fa2 theta2] = ArrayFactor_LinearPhaseDist(N,a,d,f, delta);
 Fa2 = Fa2 ./max(Fa2); % normalize the array

plot(theta1, 10*log10(Fa1), theta2, 10*log10(Fa2))
xlabel('Angle \theta (deg)')
ylabel('Array Factor (dB)')
legend('Constant Amplitude', 'Binomial Distribution',3)
grid
axis([ -180 180 -40 0])
    