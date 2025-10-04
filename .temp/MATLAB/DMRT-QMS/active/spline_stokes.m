function [Sd,Su,Sd0,Su0]=spline_stokes(nl,mu,Kprime,k)
% mu = -mu(1:ndeg);
ndeg = length(mu);

epsilon = (Kprime/k) .^2;                  %permittivity of dry snow
Sd = zeros(8*ndeg,8*ndeg,nl-1);            %Sd:S12,S23,...Sn-1,n Su:S21,S32,...Sn,n-1
Su = zeros(8*ndeg,8*ndeg,nl-1);
Sd0 = zeros(4*ndeg,4*ndeg,nl-1);
Su0 = zeros(4*ndeg,4*ndeg,nl-1);

for ii=1:nl-1
    vd = spline16(epsilon(ii),epsilon(ii+1),mu);
    vu = spline16(epsilon(ii+1),epsilon(ii),mu);
    for lin=0:ndeg - 1
        for col=0:ndeg - 1
            for jj=1:8
                Sd(8*lin+jj,8*col+jj,ii) = vd (lin+1,col+1);
                Su(8*lin+jj,8*col+jj,ii) = vu (lin+1,col+1);
            end
            for kk=1:4
                Sd0(4*lin+kk,4*col+kk,ii) = vd (lin+1,col+1);
                Su0(4*lin+kk,4*col+kk,ii) = vu (lin+1,col+1);
            end
        end
    end
end

