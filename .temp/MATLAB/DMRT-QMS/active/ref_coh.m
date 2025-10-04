function [rng0, rng] = ref_coh(rc,tai,tagd)
% interpolate the coherent reflection coefficient into gridded thetai, and
% construct the reflection matrix as required in the rough boundary
% condition.

[~,~,ntai] = size(rc);
ntagd = length(tagd);

% interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RC = zeros(4,4,ntagd);
for ii = 1:4
    for jj = 1:4
        RC(ii,jj,:) = spline(tai,rc(ii,jj,:),tagd);
        if sn(ii,jj) == ntai % positive
            RC(ii,jj,:) = min(max(RC(ii,jj,:),0),1);
        elseif sn(ii,jj) == 0 % negative or zero
            RC(ii,jj,:) = max(min(RC(ii,jj,:),0),-1);
        end
    end
end

rng0 = zeros(4*ntagd,4*ntagd);
rng = zeros(8*ntagd,8*ntagd);
for ideg = 1:ntagd
    id = block(ideg,4);
    rng0(id,id) = RC(:,:,ideg);
    id = block(ideg,8);
    rng(id(1:4),id(1:4)) = RC(:,:,ideg);
    rng(id(5:8),id(5:8)) = RC(:,:,ideg);
end

end

function id = block(x,s)
id = (x-1)*s + 1:x*s;
end