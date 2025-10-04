function [Au,Ad,B,A1,An,B1,Bn] = FillBC(d,eigvp,eigvm,Au,Ad,A1,An,...
    ke,ctinc,part1d,part1u,part2d,part2u,...
    rdown,rup,tdown,tup,Sd,Su)
% Fill the element for the boundary condition problem.

nl = length(d);
[dimh,~] = size(eigvp) ;
dimm = 2*dimh;

rS = 1:dimh;
rB = dimh + 1:dimm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill Au, Ad and B.
% Au=zeros(dimm,dimm,nl-1);
% Ad=zeros(dimm,dimm,nl-1);
BB=zeros(dimm,2,nl-1);
BC=zeros(dimm,2,nl-1);
CC=zeros(dimm,2,nl-1);
% B=zeros(dimm,2,nl-1);

for nk=1:nl-1  %1999
    % attenuate Ad and Au
    for ia = 1:dimh
        arg = eigvp(ia,nk)*d(nk);
        if(abs(real(arg)) < 20)
            Ad(:,ia,nk) = Ad(:,ia,nk) * exp(-arg);
        else
            Ad(:,ia,nk) = 0;
        end
        arg = eigvm(ia,nk+1)*d(nk+1);
        if(abs(real(arg)) < 20)
            Au(:,ia + dimh,nk) = Au(:,ia + dimh,nk)*exp(arg);
        else
            Au(:,ia + dimh,nk) = 0;
        end
    end

    BB(rS,:,nk) = part1d(:,:,nk);
    BB(rB,:,nk) = part1d(:,:,nk);
    BC(rS,:,nk) = -part1u(:,:,nk);
    arg = ke(nk) * d(nk) / ctinc(nk);
    if (abs(real(arg))<20)
        BB(rS,:,nk) = BB(rS,:,nk) + part2d(:,:,nk) * exp(-arg);
        BB(rB,:,nk) = BB(rB,:,nk) + part2d(:,:,nk) * exp(-arg);
        BC(rS,:,nk) = BC(rS,:,nk) - part2u(:,:,nk) * exp(-arg);
    end
    BB(rS,:,nk) = rdown(:,:,nk) * BB(rS,:,nk);
    BB(rB,:,nk) = -Sd(:,:,nk) * tdown(:,:,nk) * BB(rB,:,nk) ;

    CC(rS,:,nk) = part2u( :, :,nk + 1);
    CC(rB,:,nk) = part2u( :, :,nk + 1);
    BC(rB,:,nk) = part2d( :, :,nk + 1);
    arg = ke(nk+1) * d(nk+1) / ctinc(nk+1);
    if (abs(real(arg))<20)
        CC(rS,:,nk) = CC(rS,:,nk) + part1u(:,:,nk+1) * exp(-arg);
        CC(rB,:,nk) = CC(rB,:,nk) + part1u(:,:,nk + 1) * exp(-arg);
        BC(rB,:,nk) = BC(rB,:,nk) + part1d(:,:,nk + 1) * exp(-arg);
    end
    CC(rS,:,nk) = Su(:,:,nk) * tup(:,:,nk + 1) * CC( rS, :,nk);
    CC(rB,:,nk) = - rup(:,:,nk + 1) * CC(rB,:,nk);
end  %1999
B = BB + CC + BC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fill A1, An, B1 and Bn.
% A1 = zeros(8*ndeg,16*ndeg);
% An = zeros(8*ndeg,16*ndeg);
% B1 = zeros(dimh,2);
% Bn = zeros(dimh,2);

% attenuation A1 and An.
for ia=1:dimh
    arg = eigvm(ia,1) * d(1);
    if( abs(real(arg))<20)
        A1( :,ia + dimh ) = A1( :,ia + dimh ) * exp(arg);
    else
        A1( :,ia + dimh ) = 0;
    end
    arg = eigvp(ia,nl) * d(nl);
    if( abs(real(arg))<20)
        An( :, ia ) = An( :, ia ) * exp(-arg);
    else
        An( :, ia ) = 0;
    end
end

B1 = rup(:,:,1) * part2u(:,:,1)-part2d(:,:,1);
arg = ke(1) * d(1) / ctinc(1);
if (abs(real(arg))<20)
    B1 = B1+( rup(:,:,1) * part1u(:,:,1)-part1d(:,:,1))*exp(-arg);
end
Bn = rdown(:,:,nl) * part1d(:,:,nl)-part1u(:,:,nl);
arg = ke(nl) * d(nl) / ctinc(nl);
if (abs(real(arg))<20)
    Bn = Bn + ( rdown(:,:,nl) * part2d(:,:,nl) - part2u(:,:,nl) ) * exp(-arg);
end

end