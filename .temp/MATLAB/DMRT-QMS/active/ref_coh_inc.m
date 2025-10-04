function [RC] = ref_coh_inc(rc,tai,tainc)
% interpolate the coherent reflection coefficient into incident thetai, all
% angles are in radian unit.
% tainc is a scalor.

[~,~,ntai] = size(rc);
% ntagd = length(tainc);

% interpolation
% check sign of each element
sn = zeros(4,4); 
for kk = 1:ntai
    for ii = 1:4
        for jj = 1:4
            if rc(ii,jj,kk) > 0
                sn(ii,jj) = sn(ii,jj) + 1;
            end
        end
    end
end

RC = zeros(4,4);
for ii = 1: 4
    for jj = 1:4
        RC(ii,jj) = spline(tai,rc(ii,jj,:),tainc);
        if sn(ii,jj) == ntai % positive
            RC(ii,jj) = min(max(RC(ii,jj),0),1);
        elseif sn(ii,jj) == 0 % negative or zero
            RC(ii,jj) = max(min(RC(ii,jj),0),-1);
        end
    end
end

end