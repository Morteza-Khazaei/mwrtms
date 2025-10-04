function coeam = SolveBC(Au,Ad,B,A1,An,B1,Bn)
% solve the boundary problem
% coeam dimension (dimm*nl,2), second dimension for v h pol.

[dimm,~,nl] = size(Au);
nl = nl + 1;
dimh = dimm/2;

rS = 1:dimh;
rB = dimh + 1:dimm;

% AAA*[a1,a2,...an]=BBB
% AAAm=zeros(dimm*nl,dimm*nl);
% BBBm=zeros(dimm*nl,2);
AAAm=sparse(dimm*nl,dimm*nl);
BBBm=sparse(dimm*nl,2);

dr = block(1,dimm);
AAAm(dr(rS),dr)=A1;
BBBm(dr(rS),:) = B1;

for ii=1:nl-1
    dr = block(ii,dimm);
    AAAm(dr + dimh,dr) = Ad(:,:,ii);
    AAAm(dr + dimh,dr + dimm) = - Au(:,:,ii);
    BBBm(dr + dimh,:) = B(:,:,ii);
end
dr = block(nl,dimm);
AAAm(dr(rB), dr)=An;
BBBm(dr(rB),:) = Bn;

coeam=AAAm\BBBm;
end


function id = block(x,s)
id = (x-1)*s + 1:x*s;
end