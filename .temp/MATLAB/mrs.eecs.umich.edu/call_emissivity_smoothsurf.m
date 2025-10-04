%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

theta = 0:1:90;

np = length(theta);

%---------------------------------------
 eps = 3.5 -1i * 0.2;

for ii = 1:np
    
    [e_v1(ii) e_h1(ii)]= Emissivity_SmoothSurf(eps, theta(ii));
end

%----------------------------------------

eps = 17.9 -1i * 7.2;
for ii = 1:np
    
    [e_v2(ii) e_h2(ii)]= Emissivity_SmoothSurf(eps, theta(ii));
end

%-----------------------------------------
 eps = 54.4 -1i * 36.8;
for ii = 1:np
    
    [e_v3(ii) e_h3(ii)]= Emissivity_SmoothSurf(eps, theta(ii));
end


plot(theta, e_h1, theta, e_v1, theta, e_h2, theta, e_v2, ...
    theta, e_h3, theta, e_v3)
axis([0 90 0 1])

xlabel('\theta_1 (deg)')
ylabel('Emissivity')
grid