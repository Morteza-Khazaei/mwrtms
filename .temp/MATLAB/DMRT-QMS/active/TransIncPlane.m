function [sigma,gamma,deg0int] = TransIncPlane(I1U0,I1U,mu,n10,tai0inc)
% Transfrom the specific intensity spectrum into intensity in the 
% incident plane (phi = pi, theta[0, pi/2]) (phi = 0, theta [0,pi/2])
% by summing up Fourier Series
% Relate the angle in region 1 and zero by Snell's law, considering
% critical angle effect.
% Also calculate backscattering by interpolation.
% 
% mu - upward angles correlating to I1U
% n10 : n1/n0 = k1/k0.
% tai0inc: incident angle
% 
% sigma, - backscattering in dB scale, (4,2), Obs: (Iv,Ih,U,V), Inc: (v,h)
% deg0int, gamma, bistatic scattering in dB scale, (4,2,ndeg0int).

ndeg = length(mu);
[~,npol,Mmax] = size(I1U);

% remove the zero imaginary part, perform calculation in real domain.
I1U0 = real(I1U0);
I1U = real(I1U);

%m=0, zeroth harmonic
Iupf = reshape(I1U0,4,ndeg,npol);    % phi = 0; dimension (4*ndeg,2)
Iupb = reshape(I1U0,4,ndeg,npol);    % phi = pi;

I1U = reshape(I1U,8,ndeg,npol,Mmax);
for m = 1:Mmax
    Iupf = Iupf + I1U(1:4,:,:,m);
    Iupb = Iupb + (-1)^m * I1U(1:4,:,:,m);
end

crit = asin(1/n10);
range = mu > cos(crit);
tai = acos(mu(range)); % valid angle in region 1
tai0 = asin(n10*sin(tai));     % valid angle in region 0
deg0 = tai0/pi*180;
deg0inc = tai0inc/pi*180;
fac = 4 * pi * cos(tai0) / cos(tai0inc); % coefficient relating I to gamma

deg0int = [- flipud(deg0(:));deg0(:)];
ndeg0 = sum(range);
gamma = zeros(4,npol,2*ndeg0);
% sigma = zeros(4,npol);
for ideg = 1:ndeg0
    Iupf(:,ideg,:) = Iupf(:,ideg,:)*fac(ideg);
    Iupb(:,ideg,:) = Iupb(:,ideg,:)*fac(ideg);
end
for ipol = 1:npol
    gamma(:,ipol,ndeg0+1:2*ndeg0) = reshape(Iupf(:,1:ndeg0,ipol),4,1,ndeg0);
    gamma(:,ipol,ndeg0:-1:1) = reshape(Iupb(:,1:ndeg0,ipol),4,1,ndeg0);
end
% gamma = gamma*fac;
sigma = spline(deg0int,gamma,-deg0inc)*cos(tai0inc);

gamma = 10*log10(gamma);
sigma = 10*log10(sigma);

end
