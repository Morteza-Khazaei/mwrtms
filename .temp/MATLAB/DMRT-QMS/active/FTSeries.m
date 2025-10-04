function [g0,gc,gs] = FTSeries(phs,integ,Mmax)

gc = zeros(Mmax,1);
gs = zeros(Mmax,1);
g0 = trapz(phs,integ)/(2*pi);
% m-th harmonics
for m = 1:Mmax
    gc(m) = trapz(phs,integ.*cos(m*phs))/pi;
    gs(m) = trapz(phs,integ.*sin(m*phs))/pi;
end

end