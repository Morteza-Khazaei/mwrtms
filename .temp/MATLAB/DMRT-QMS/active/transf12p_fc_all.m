function [P0,PC,PS] = transf12p_fc_all(bmua,P11,P22,P12,P21,P33,P34,mui,mus,nphs,Mmax,...
                        P13,P14,P23,P24,P31,P32,P41,P42)
% Function transforms the Phase Matrix P from 1-2 frame to main frame, and
% then perform Fourier Transform with respect to azimuth angle phi. Phase
% matrix P is defined on angle bmua. The main frame space is sampled by mui
% and mus. nphs is the number of azimuth angles in the main frame to
% facilate the Fourier Transform.
% Default value: nphs = 61; Mmax = 4; 
% off diagal phase matrix elements P13,P14,P23,P24,P31,P32,P41,P42 are
% assumed to be zero if not given.
% 

np = length(bmua);
ntai = length(mui);
ntas = length(mus);

if nargin < 12
   P13 = zeros(1,np); 
   P14 = zeros(1,np); 
   P23 = zeros(1,np); 
   P24 = zeros(1,np); 
   P31 = zeros(1,np); 
   P32 = zeros(1,np); 
   P41 = zeros(1,np); 
   P42 = zeros(1,np); 
   if nargin < 11
       Mmax = 4;
       if nargin < 10
           nphs = 61;
       end
   end
end

tai = acos(mui);
cti = mui;
sti = sin(tai);
tas = acos(mus);
cts = mus;
sts = sin(tas);
phs=linspace(0,2*pi,nphs);
cps = cos(phs);
sps = sin(phs);

% relate direction in main frame with \Theta in 1-2 frame
npt = ntai*ntas*nphs;
cbt = zeros(1,npt);
for idegi = 1:ntai
    n0i = (idegi - 1)*ntas*nphs;
    for idegs = 1:ntas
        n0s = n0i + (idegs - 1)*nphs;
        cbt((1:nphs) + n0s) = cti(idegi)*cts(idegs) + sti(idegi)*sts(idegs)*cps;
    end
end

% interpolate
P11i = interp1(bmua,P11,cbt,'spline');
P22i = interp1(bmua,P22,cbt,'spline');
P12i = interp1(bmua,P12,cbt,'spline');
P21i = interp1(bmua,P21,cbt,'spline');
P33i = interp1(bmua,P33,cbt,'spline');
P34i = interp1(bmua,P34,cbt,'spline');

P13i = interp1(bmua,P13,cbt,'spline');
P14i = interp1(bmua,P14,cbt,'spline');
P23i = interp1(bmua,P23,cbt,'spline');
P24i = interp1(bmua,P24,cbt,'spline');
P31i = interp1(bmua,P31,cbt,'spline');
P32i = interp1(bmua,P32,cbt,'spline');
P41i = interp1(bmua,P41,cbt,'spline');
P42i = interp1(bmua,P42,cbt,'spline');

% transform frame
PA = zeros(ntas,nphs,ntai,4,4);

hi = [0 1 0];
for idegi = 1:ntai
    n0i = (idegi - 1)*ntas*nphs;
    ki = [sti(idegi) 0 cti(idegi)];
    vi = [cti(idegi) 0 -sti(idegi)];
    for idegs = 1:ntas
        n0s = n0i + (idegs - 1)*nphs;
        for iphs = 1:nphs
            n0 = n0s + iphs;
            
            ks = [sts(idegs)*cps(iphs) sts(idegs)*sps(iphs) cts(idegs)];
            vs = [cts(idegs)*cps(iphs) cts(idegs)*sps(iphs) -sts(idegs)];
            hs = [-sps(iphs) cps(iphs) 0];
            
            onei = [0 -1 0]; % forward or backward direction.
            q = cross3(ks,ki);
            absq2 = dot3(q,q);
            if absq2 > 1e-9
                onei = q/sqrt(absq2);
            end
            oness = onei;
            % twoi=cross(ki,onei);
            % twos=cross(ks,oness);
            alphai = atan2(dot3(hi,onei),dot3(vi,onei));
            alphas = atan2(dot3(hs,oness),dot3(vs,oness ));
            s2i = sin(alphai)^2;
            c2i = 1 - s2i; % cos(alphai)^2;
            s2ai = sin(2*alphai);
            c2ai = c2i - s2i; % cos(2*alphai);
            s2s = sin(alphas)^2;
            c2s = 1 - s2s; % cos(alphas)^2;
            s2as=sin(2*alphas);
            c2as= c2s - s2s; % cos(2*alphas);
            Pti=[c2i s2i s2ai/2 0;s2i c2i -s2ai/2 0;-s2ai s2ai c2ai 0;0 0 0 1];
            Pts=[c2s s2s -s2as/2 0;s2s c2s s2as/2 0;s2as -s2as c2as 0;0 0 0 1];
            
            P12f=[  P11i(n0) P12i(n0) P13i(n0) P14i(n0);
                    P21i(n0) P22i(n0) P23i(n0) P24i(n0);
                    P31i(n0) P32i(n0) P33i(n0) P34i(n0);
                    P41i(n0) P42i(n0) -P34i(n0) P33i(n0); ];
            PA(idegs,iphs,idegi,:,:) = Pts*P12f*Pti;
        end
    end
end


% Fourier expansions
P0=zeros(ntas,ntai,4,4);
PC=zeros(Mmax,ntas,ntai,4,4);
PS=zeros(Mmax,ntas,ntai,4,4);

for idegi = 1:ntai  %200
    for idegs = 1:ntas  %201
        % zero-th harmonic
        for ii = 1:4
            for jj = 1:4
                integ = reshape(PA(idegs,:,idegi,ii,jj),1,nphs)/(2*pi);
                P0(idegs,idegi,ii,jj)=trapz(phs,integ);
            end
        end
        % m-th harmonic
        for m = 1:Mmax  %202
            cphm = cos(m*phs)/pi;
            sphm = sin(m*phs)/pi;
            for ii = 1:4  %203
                for jj = 1:4 %204
                    integc = reshape(PA(idegs,:,idegi,ii,jj),1,nphs).*cphm;
                    integs = reshape(PA(idegs,:,idegi,ii,jj),1,nphs).*sphm;
                    PC(m,idegs,idegi,ii,jj)=trapz(phs,integc);
                    PS(m,idegs,idegi,ii,jj)=trapz(phs,integs);
                end %204
            end %203
        end %202
    end %201
end %200

end

function c = cross3(a,b)
c1 = a(2)*b(3) - a(3)*b(2);
c2 = a(3)*b(1) - a(1)*b(3);
c3 = a(1)*b(2) - a(2)*b(1);
c = [c1 c2 c3];
end

function c = dot3(a,b)
c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
end