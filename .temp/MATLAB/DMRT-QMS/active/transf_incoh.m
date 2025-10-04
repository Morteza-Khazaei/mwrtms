function [gammaic0t, gammaic0tinc, gammaict, gammaictinc] = ...
    transf_incoh(gammainc,tas,phs,tai,tagd,taiinc,Mmax)
% transform the incoherent reflection coefficient into Fourier series with 
% respect to phs, and resample at the grid point of theta. 
% Finally, put variables in order as required by the boundary conditions
% all angles in radian
% phi = 0, incident phi fixed to zero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % approach one, first carry out Fourier transform, 
% % then resample to grid points in Fourier Domain.
% [~,~,ntas,nphs,ntai] = size(gammainc);
% G0 = zeros(4,4,ntas,ntai);
% GC = zeros(4,4,ntas,ntai,Mmax);
% GS = zeros(4,4,ntas,ntai,Mmax);
% 
% % carry out Fourier transform
% for idegi = 1:ntai
%     for idegs = 1:ntas
%         for ii = 1:4
%             for jj = 1:4
%                 % zero-th harmonic
%                 integ = reshape(gammainc(ii,jj,idegs,:,idegi),1,nphs);
% %                 G0(ii,jj,idegs,idegi) = trapz(phs,integ)/(2*pi);
% %                 % m-th harmonics
% %                 for m = 1:Mmax
% %                     GC(ii,jj,idegs,idegi,m) = trapz(phs,integ.*cos(m*phs))/pi;
% %                     GS(ii,jj,idegs,idegi,m) = trapz(phs,integ.*sin(m*phs))/pi;
% %                 end
%                 [G0(ii,jj,idegs,idegi),GC(ii,jj,idegs,idegi,:),GS(ii,jj,idegs,idegi,:)] ...
%                     = FTSeries(phs,integ,Mmax);
%             end
%         end
%     end
% end
% 
% % resample to grid points of (thetas, thetai) and (thetainc,thetai)
% ntagd = length(tagd);
% G0GD = zeros(4,4,ntagd,ntagd);      % thetas by thetai
% GCGD = zeros(4,4,ntagd,ntagd,Mmax);
% GSGD = zeros(4,4,ntagd,ntagd,Mmax);
% 
% G0inc = zeros(4,4,ntagd);
% GCinc = zeros(4,4,ntagd,Mmax);
% GSinc = zeros(4,4,ntagd,Mmax);
% 
% [XI,YI] = meshgrid(tagd,tagd); % (thetai,thetas)
% [xi,yi] = meshgrid(tagd,taiinc); % (thetai, thetas = thetainc)
% for ii = 1:4
%     for jj = 1:4
%         % zero-th harmonic
%         G0GD(ii,jj,:,:) = interp2(tai,tas,reshape(G0(ii,jj,:,:),ntas,ntai),XI,YI,'*spline');
%         G0inc(ii,jj,:) = interp2(tai,tas,reshape(G0(ii,jj,:,:),ntas,ntai),xi,yi,'*spline');
%         % m-th harmonic
%         for m = 1:Mmax
%             GCGD(ii,jj,:,:,m) = interp2(tai,tas,reshape(GC(ii,jj,:,:,m),ntas,ntai),XI,YI,'*spline');
%             GSGD(ii,jj,:,:,m) = interp2(tai,tas,reshape(GS(ii,jj,:,:,m),ntas,ntai),XI,YI,'*spline');
%             GCinc(ii,jj,:,m) = interp2(tai,tas,reshape(GC(ii,jj,:,:,m),ntas,ntai),xi,yi,'*spline');
%             GSinc(ii,jj,:,m) = interp2(tai,tas,reshape(GS(ii,jj,:,:,m),ntas,ntai),xi,yi,'*spline');
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second approach, first resampe to grid points wrt. tas and tai,
% then carry out Fourier transform.
[~,~,ntas,nphs,ntai] = size(gammainc);
ntagd = length(tagd);

% interpolation
G = zeros(4,4,ntagd,nphs,ntagd);
Ginc = zeros(4,4,nphs,ntagd);

[XI,YI] = meshgrid(tagd,tagd); % (thetai,thetas)
% [xi,yi] = meshgrid(tagd,taiinc); % (thetai, thetas = thetainc)
[xi,yi] = meshgrid(taiinc,tagd); % (thetai = thetainc, thetas), resiprocity

for kk = 1:nphs
    for ii = 1:4 % 4
        for jj = 1:4 %4
            ip = interp2(tai,tas,reshape(gammainc(ii,jj,:,kk,:),ntas,ntai),XI,YI,'*spline');
            ipi = interp2(tai,tas,reshape(gammainc(ii,jj,:,kk,:),ntas,ntai),xi,yi,'*spline');
            if ii < 3 && jj < 3
                ip = min(max(ip,0),1);
                ipi = min(max(ipi,0),1);
            end
            G(ii,jj,:,kk,:) = reshape(ip,1,1,ntagd,1,ntagd);
            Ginc(ii,jj,kk,:) = reshape(ipi,1,1,1,ntagd);
        end
    end
end

% Fourier Transform
G0GD = zeros(4,4,ntagd,ntagd);      % thetas by thetai
GCGD = zeros(4,4,ntagd,ntagd,Mmax);
GSGD = zeros(4,4,ntagd,ntagd,Mmax);

G0inc = zeros(4,4,ntagd);
GCinc = zeros(4,4,ntagd,Mmax);
GSinc = zeros(4,4,ntagd,Mmax);

for idegi = 1:ntagd
    for idegs = 1:ntagd
        for ii = 1:4
            for jj = 1:4
                integ = reshape(G(ii,jj,idegs,:,idegi),1,nphs);
                [G0GD(ii,jj,idegs,idegi),GCGD(ii,jj,idegs,idegi,:),GSGD(ii,jj,idegs,idegi,:)] ...
                    = FTSeries(phs,integ,Mmax);
            end
        end
    end
end

for idegi = 1:ntagd % idegs for resiprocity
    for ii = 1:4
        for jj = 1:4
           integ_inc = reshape(Ginc(ii,jj,:,idegi),1,nphs);
            [G0inc(ii,jj,idegi),GCinc(ii,jj,idegi,:),GSinc(ii,jj,idegi,:)] ...
                = FTSeries(phs,integ_inc,Mmax);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% place in order as required by the boundary conditions
gammaic0t = zeros(4*ntagd,4*ntagd);
gammaic0tinc = zeros(4*ntagd,4);
gammaict = zeros(8*ntagd,8*ntagd,Mmax);
gammaictinc = zeros(8*ntagd,4,Mmax);
for idegi = 1:ntagd
    for idegs = 1:ntagd
        idx = block(idegi,4);
        idy = block(idegs,4);
        gammaic0t(idx,idy) = G0GD(:,:,idegs,idegi)';
        idx = block(idegi,8);
        idy = block(idegs,8);
        for m = 1:Mmax
            gammaict(idx(1:4),idy(1:4),m) = GCGD(:,:,idegs,idegi,m)';
            gammaict(idx(1:4),idy(5:8),m) = GSGD(:,:,idegs,idegi,m)';
            gammaict(idx(5:8),idy(1:4),m) = - GSGD(:,:,idegs,idegi,m)';
            gammaict(idx(5:8),idy(5:8),m) = GCGD(:,:,idegs,idegi,m)';
        end
    end
end
for idegi = 1:ntagd
    coe = cos(taiinc)/cos(tagd(idegi));
    idx = block(idegi,4);
%     gammaic0tinc(idx,1:4) = G0inc(:,:,idegi)';
    gammaic0tinc(idx,1:4) = G0inc(:,:,idegi)*coe; %resiprocity
    idx = block(idegi,8);
    for m = 1:Mmax
%         gammaictinc(idx(1:4),1:4) = GCinc(:,:,idegi,m)';
%         gammaictinc(idx(5:8),1:4) = - GSinc(:,:,idegi,m)';
        gammaictinc(idx(1:4),1:4) = GCinc(:,:,idegi,m)*coe; %resiprocity
        gammaictinc(idx(5:8),1:4) = GSinc(:,:,idegi,m)*coe; %resiprocity
    end
end

end

function id = block(x,s)
id = (x-1)*s + 1:x*s;
end