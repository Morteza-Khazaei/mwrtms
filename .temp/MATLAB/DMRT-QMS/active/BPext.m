function [A1, An] = BPext(vecupt,vecdpt,vecumt,vecdmt,rupt,...
    vecupb,vecdpb,vecumb,vecdmb,rdownb)
% preliminary job for the boundary condition problem.
% external boundaries (top and bottom)

% ndeg
[dimh,~] = size(vecupt); % 8*ndeg
dimm = 2*dimh;
rS = 1:dimh;
rB = dimh + 1:dimm;

A1=zeros(dimh,dimm);
An=zeros(dimh,dimm);

% fill out matrix for the simulataneous equation. 
% upward
A1(rS,rS)  = vecdpt - rupt*vecupt;
A1(rS,rB)  = vecdmt - rupt*vecumt;    % to be att by exp(eigvm(ia,1) * d(1))
% downward
An(rS,rS)  = vecupb - rdownb*vecdpb;    % to be att by exp(-eigvp(ia,nl) * d(nl))
An(rS,rB)  = vecumb - rdownb*vecdmb;

end
