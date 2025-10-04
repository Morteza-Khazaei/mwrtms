function Bs = BCrs(ke,d,ctinc,part1d,part2d,N_i,mu,aa,gammaict,gammaictinc)
% modify Bn due to rough surface effect
% ke,d,ctinc,part1d,part2d and N_i, are of the last layer.

[dimm,~] = size(gammaict); %8*ntagd
ndeg2 = length(mu); %ndeg2
ndeg = ndeg2/2;
mut = mu(1:ndeg);
Nel = dimm/ndeg2; %8

dimh = dimm/2;

Bs0 = zeros(dimh,2);
Bs1 = zeros(dimh,2);
Bs2 = zeros(dimh,2);

for ipol = 1:2
    % modify Bn
    % vec = weight(part1d(:,nl),aa(1:ndeg),8);
    % Bs1 = gammaict(:,:,m)*vec/4;
    Bs1(:,ipol) = weight_int(part1d(:,ipol),mut,gammaict(:,:),(mu + 1)/2,aa,Nel)/(Nel/2);
    % Bn(1:8*ndeg) = Bn(1:8*ndeg) + Bs1;
    arg = ke*d/ctinc;
    if(abs(real(arg)) < 20)
      % vec = weight(part2d(:,nl),aa(1:ndeg),8);
      % Bs2 = gammaict(:,:,m)*vec/4*exp(-arg);
      %Bs0 = gammaictinc(:,:,m)*N_i(:,nl,ipol)/4/pi*exp(-arg);
      Bs2(:,ipol) = weight_int(part2d(:,ipol),mut,gammaict(:,:),(mu + 1)/2,aa,Nel)/(Nel/2)*exp(-arg);
      Bs0(:,ipol) = multiply_folder(N_i(:,ipol),mut,gammaictinc(:,:),(mu + 1)/2,Nel)/4/pi*exp(-arg);
%       Bn(1:8*ndeg) = Bn(1:8*ndeg) + Bs0 + Bs2;
    end
end
Bs = Bs0 + Bs1 + Bs2;

end