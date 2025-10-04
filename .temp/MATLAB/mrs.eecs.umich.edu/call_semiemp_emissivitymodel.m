%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

s = 0.01;
l = 0.1;

theta = 0:1:80;
np = length(theta);

eps = (9 - 1i*0.9);

for n = 1:np
    [e_v(n) e_h(n)] = SemiEmp_EmissivityModel_RoughSurf(eps, s,l,theta(n));
end

figure(1)

plot(theta, e_v, '-r', theta, e_h,'-b')
xlabel('\theta (deg)')
ylabel('Emissivity')
legend('e_v', 'e_h', 3)
grid
