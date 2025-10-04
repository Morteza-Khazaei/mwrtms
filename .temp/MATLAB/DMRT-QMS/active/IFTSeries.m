function [integ] = IFTSeries(phs,g0,gc,gs)

nphs = length(phs);
Mmax = length(gc);
integ = zeros(nphs,1);
for i = 1:nphs
    mphs = phs(i)*(1:Mmax)';
    cs = cos(mphs);
    ss = sin(mphs);
    integ(i) = g0 + sum(gc.*cs + gs.*ss);
end

end