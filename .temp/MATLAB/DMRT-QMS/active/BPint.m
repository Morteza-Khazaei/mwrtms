function [Au, Ad] = BPint(vecup,vecdp,vecum,vecdm,rdown,tdown,Sd,...
    vecupn,vecdpn,vecumn,vecdmn,rupn,tupn,Su)
% preliminary job for the boundary condition problem.
% internal boundaries

% ndeg
[dimh,~] = size(vecup); % 8*ndeg
dimm = 2*dimh;
rS = 1:dimh;
rB = dimh + 1:dimm;

Au=zeros(dimm,dimm);
Ad=zeros(dimm,dimm);

% fill out matrix for the simulataneous equation. 
% downward
Ad(rS,rS)  = vecup - rdown*vecdp;  % to be att by exp(-eigvp(ia,nk)*d(nk))
Ad(rB,rS)  = Sd*tdown*vecdp;       % to be att
Ad(rS,rB)  = vecum - rdown*vecdm;
Ad(rB,rB)  = Sd*tdown*vecdm;

% upward
Au(rS,rS)  = Su*tupn*vecupn;
Au(rB,rS)  = vecdpn - rupn*vecupn;
Au(rS,rB)  = Su*tupn*vecumn;       % to be att by exp(eigvm(ia,nk+1)*d(nk+1))
Au(rB,rB)  = vecdmn - rupn*vecumn; % to be att

end
