function [As] = BPrs(gammaict,mu,aa,vecdmb,vecdpb)
% modify An and Bn, due to rough snow-ground surface

[dimm,~] = size(gammaict); %8*ntagd
ndeg2 = length(mu); %ndeg2
ndeg = ndeg2/2;
mut = mu(1:ndeg);
Nel = dimm/ndeg2; %8

dimh = dimm/2;

As = zeros(dimh,dimm);
for ia = 1:dimh
    % vec = weight(vecdp(:,ia,nl),aa(1:ndeg),8);
    % Asp = gammaict(:,:,m)*vec/4*exp(-arg);
    As(:,ia)        = weight_int(vecdpb(:,ia),mut,gammaict,(mu + 1)/2,aa,Nel)/(Nel/2); % to be att by exp(-eigvp(ia,nl)*d(nl))
    
    % vec = weight(vecdm(:,ia,nl),aa(1:ndeg),8);
    % Asm = gammaict(:,:,m)*vec/4;
    As(:,ia + dimh) = weight_int(vecdmb(:,ia),mut,gammaict,(mu + 1)/2,aa,Nel)/(Nel/2);
end

end