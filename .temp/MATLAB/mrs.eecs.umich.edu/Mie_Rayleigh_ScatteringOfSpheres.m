%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.12: Mie and Rayleigh Scattering By Spherical Particle

%Description: Code computes the absorption, scattering, extinction, and
%backscattering efficiencies of a dielectric sphere according to bot the
%Mie solution and the Rayleigh approximation.

%Input Variables
    %r: radius of particle (meters)
    %f: frequency (GHz)
    %epsp: epsp1-j*epsp2 for particles
    %epsb = epsb1 -j epsb2 for background medium 
    
%Output Variables
    %Es: Mie Scattering Efficiency 
    %Ea: Mie Absorption Efficiency
    %Ee: Mie Extinction Efficiency
    %Eb: Mie Backscattering Efficiency
    
    %eta_s_r: Rayleigh Scattering Efficiency 
    %eta_a_r: Rayleigh Absorption Efficiency
    %eta_e_r: Rayleigh Extinction Efficiency
    %eta_b_r: Rayleigh Backscattering Efficiency

%Matlab Code

function [Es Ea Ee Eb eta_s_r eta_a_r eta_e_r eta_b_r] = ... 
    Mie_Rayleigh_ScatteringOfSpheres(r, f, epsp, epsb)

epsb1 = real(epsb);

np = sqrt(epsp); % index of refraction of spherical particle
nb = sqrt(epsb); % index of refraction of background medium

n = np ./nb; % relative index of refraction
%  n = np ./sqrt(epsb1);


chi = 20./3 *pi*r*f*sqrt(epsb1); % normalized circumference in background

%--- Calculate Rayleigh Approximation solution

BigK = (n.^2 - 1) ./(n.^2 + 2); 

eta_s_r = 8/3 .* chi.^4 .* abs(BigK).^2;

eta_a_r = 4 .* chi .* imag(-BigK); 

eta_e_r = eta_a_r + eta_s_r; % Extinction efficiency

eta_b_r = 4 * chi.^4 .* abs(BigK).^2;  % backscattering efficiency



%--- Calculate Mie Scattering solution

%Calculation of Es
l=1;
first = true;
runSum=linspace(0,0,numel(f));
oldSum=linspace(0,0,numel(f));

    %Values of W0 and W-1
    W_1 = sin(chi)+1i*cos(chi);
    W_2 = cos(chi)-1i*sin(chi);
    
    %Value of A0
    A_1 = cot(n*chi);
            
    %A_1=    (sin(real(n).*chi).*cos(real(n).*chi)+j*sinh(imag(n).*chi).*cosh(imag(n).*chi))...
    %        ./(sin(real(n).*chi).^2 + sinh(imag(n).*chi).^2);

while (first || ~endSum(oldSum, runSum, numel(chi)))
    W=(2*l-1)./chi .* W_1 - W_2;

    A = -l./(n.*chi) + (l./(n.*chi)-A_1).^-1;
    
    a = ((A/n + l./chi).*real(W)-real(W_1)) ./ ((A/n+l./chi).*W-W_1);
    b = ((n*A + l./chi).*real(W)-real(W_1)) ./ ((n*A+l./chi).*W-W_1);
    
    sumTerm = (2*l + 1).*(abs(a).^2+abs(b).^2);
    oldSum = runSum;
    runSum = runSum + sumTerm;
    
    %Increment Index varible
    l=l+1;
    
    %Increment W terms
    W_2 = W_1;
    W_1 = W;
    
    %Increment A Terms
    A_1=A;
    
    %Set first pass to false
    first = false;
end
Es = 2./(chi).^2 .* runSum;

%Calculation of Ee
l=1;
first = true;
runSum=linspace(0,0,numel(f));
oldSum=linspace(0,0,numel(f));

    %Values of W0 and W-1
    W_1 = sin(chi)+1i*cos(chi);
    W_2 = cos(chi)-1i*sin(chi);
    
    %Value of A0
    A_1 = cot(n*chi);
            
    %A_1 =   (sin(real(n).*chi).*cos(real(n).*chi)+j*sinh(imag(n).*chi).*cosh(imag(n).*chi))...
    %       ./(sin(real(n).*chi).^2 + sinh(imag(n).*chi).^2);

while (first || ~endSum(oldSum, runSum, numel(chi)))
    W=(2*l-1)./chi .* W_1 - W_2;

    A = -l./(n.*chi) + (l./(n.*chi)-A_1).^-1;
    
    a = ((A/n + l./chi).*real(W)-real(W_1)) ./ ((A/n+l./chi).*W-W_1);
    b = ((n*A + l./chi).*real(W)-real(W_1)) ./ ((n*A+l./chi).*W-W_1);
    
    sumTerm = (2*l + 1).*real(a+b);
    oldSum = runSum;
    runSum = runSum + sumTerm;
    
    %Increment Index varible
    l=l+1;
    
    %Increment W terms
    W_2 = W_1;
    W_1 = W;
    
    %Increment A Terms
    A_1=A;
    
    %Set first pass to false
    first = false;
end
Ee = 2./(chi).^2 .* runSum;

%Calculation of Eb
l=1;
first = true;
runSum=linspace(0,0,numel(f));
oldSum=linspace(0,0,numel(f));

    %Values of W0 and W-1
    W_1 = sin(chi)+1i*cos(chi);
    W_2 = cos(chi)-1i*sin(chi);
    
    %Value of A0
    A_1 = cot(n*chi);
            

while (first || ~endSum(oldSum, runSum, numel(chi)))
    W=(2*l-1)./chi .* W_1 - W_2;

    A = -l./(n.*chi) + (l./(n.*chi)-A_1).^-1;
    
    a = ((A/n + l./chi).*real(W)-real(W_1)) ./ ((A/n+l./chi).*W-W_1);
    b = ((n*A + l./chi).*real(W)-real(W_1)) ./ ((n*A+l./chi).*W-W_1);
    
    sumTerm = (-1)^l .* (2*l + 1) .*(a-b);
    oldSum = runSum;
    runSum = runSum + sumTerm;
    
    %Increment Index varible
    l=l+1;
    
    %Increment W terms
    W_2 = W_1;
    W_1 = W;
    
    %Increment A Terms
    A_1=A;
    
    %Set first pass to false
    first = false;
end

Eb = 1./(chi).^2 .* abs(runSum).^2;


Ea=Ee-Es;

end



function [stop]=endSum(A0, A1, num)
    stop=true;
    pDiff = abs((A1-A0)./A0) .* 100;
    
    for t=1:num,
        if((pDiff(t)>=0.001) || (A0(t) ==0))
            stop=false;
        end
    end

end