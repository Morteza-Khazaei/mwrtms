function [sigma,gamma,deg0int] ...
    = AMLDRS(deg0inc, depth, perm_gnd, rsfile, ndeg, Mmax,...
    perm_eff_a,ka,ks,bmua,P11_a,P22_a,P33_a,P34_a,P12_a,P21_a,...
    P13_a, P14_a, P23_a, P24_a, P31_a, P32_a, P41_a, P42_a )
% AMLDRS (Active Multiple Scattering multi-Layer DMRT combined with Rough 
% Surface boundary condtion) calculates the backscatter and bi-static
% scatering coefficient of a planar layered snowpack laid on top of either
% a planar ground or rough ground. 
% 
% Input Arguments:
%   deg0inc - incident angle in region 0 (air) in degree, scalar.
%   depth - depth of each snow layer (a total of nl layers) from top to 
%       bottom, and multiple depth configuration could be supplied by 
%       different rows. Depth of of dimension (nconfig, nlayer). Unit in
%       cm.
%   perm_gnd - permittivity of flat ground, can be complex. If ground is
%       assumed to be rough, this para doesn't take effect, use empty [].
%   rsfile - path (string) of the rough surface charaterization file, in
%       which the coherent reflection coefficient and incoherent bi-static
%       scattering coefficients of the snow-ground interface are supplied. 
%       This para doesn't take effect if a planar ground is assume, and
%       could be leave as [].
%       rsfile is a *.mat file containing the following 5 variables:
%       'tai_rs', 'tas_rs', 'phs_rs', 'rc_rs', 'gammainc_rs'
%           tai_rs (1,ntai): incident theta angle array in radian
%           tas_rs (1,ntas): scatterig theta angle array in radian
%           phs_rs (1,nphs): scattering phi angle array in radian
%           rc_rs (4,4,ntai): coherent reflection coefficient of rs
%           gammainc_rs (4,4,ntas,nphs,ntai): incoherent bistatic
%               scattering coefficient of rough surface (rs)
%   ndeg - half of the number of quadrature angles in the discrete 
%       coordinate of theta, (i.e. the number of theta sampling in (0,pi/2))
%       default value is 16, or choose some value larger to test
%       convergence.
%   Mmax - order of Fourier harmonics expansion. Default value is 4 for
%       flat ground and 8 for rough surface, or choose some larger value to
%       test convergence.
%   perm_eff_a - effective permittivities of the snow medium in each layer,
%       should be real (the effect of the imaginary part is considered by 
%       kappaa), and be a colum vector (nl, 1). If the top halfspace is not
%       air, this should be the ratio between snow permittivites and top
%       halfspace permittivity.
%   ka,ks - absorption coefficient kappaa and scattering coefficient kappas
%       of the random dense medium (snow), should be real, and in colum
%       vectors of (nl, 1). Unit in (1/cm).
%   bmua - cos(\Theta) of phase matrix, (1, np) row vector. No restriction
%       on the length np. But np > 32 is suggested to have fine resolution
%       of phase matrix in 1-2 frame. Default value is 65.
%   Pmn - Phase matrixs NORMALIZED by kappas. The off-diagonal elements
%       P13, P14, P23, P24, P31, P32, P41, P42 are automatically filled to zere
%       if not provided. Each P is a (nl,np) matrix.
% 
% Output Parameters:
%   sigma - backscatter in dB scale, (4,nd). 
%       first dimension (colum): vv, hv, vh, hh; 
%       second dimension (row): different depth configurations (nd).
%   gamma - bistatic scattering in dB scale, (4,ndeg0,nd) 
%       first dimension: vv, hv, vh, hh;
%       second dimension: observation angles in region 0 (air) in degree,
%       in increasing order; 
%       third dimension: different layer depth
%       configurations.
% 
% Usage:
%     deg0inc = 40; % degree
%     depth = [10 20; 20 40]; % cm
%     [nd, nl] = size(depth);
%     perm_gnd = 3.2 + 0.002i;
%     rsfile = []; % assume flat ground
%     ndeg = 16;
%     Mmax = 4;
%     perm_eff_a = [1.2 1.3]';
%     ks = [0.5 1.0]'/100;    % 1/cm
%     albedo = 0.9;
%     ka = ks/albedo*(1 - albedo);
%     np = 65;
%     bmua = linspace(-1,1,np);
%     % assume simple Rayleigh phase matrix.
%     P11 = ones(1,np);
%     P22 = bmua.^2;
%     P33 = bmua;
%     ksn = pi*trapz(bmua,P11 + P22);
%     P11 = P11/ksn;
%     P22 = P22/ksn;
%     P33 = P33/ksn;
%     P11_a = zeros(nl,np);
%     P22_a = zeros(nl,np);
%     P12_a = zeros(nl,np);
%     P21_a = zeros(nl,np);
%     P33_a = zeros(nl,np);
%     P34_a = zeros(nl,np);
%     for il = 1:nl
%         P11_a(il,:) = P11;
%         P22_a(il,:) = P22;
%         P33_a(il,:) = P33;
%     end
%     [sigma,gamma,deg0int] ...
%         = AMLDRS(deg0inc, depth, perm_gnd, rsfile, ndeg, Mmax,...
%         perm_eff_a,ka,ks,bmua,P11_a,P22_a,P33_a,P34_a,P12_a,P21_a);  
%     % show the backscattering
%     for il = 1:nl
%         fprintf('Config %d: %.2f(vv) %.2f(hv) %.2f(vh) %.2f(hh)\n', ...
%             il, sigma(1,il),sigma(2,il),sigma(3,il),sigma(4,il));
%     end
%     % show the bi-static scattering for the first depth configuration.
%     figure; hold on;
%     plot(deg0int,gamma(1,:,1),'r.-', 'linewidth', 2);
%     plot(deg0int,gamma(2,:,1),'r:', 'linewidth', 2);
%     plot(deg0int,gamma(3,:,1),'b--', 'linewidth', 2);
%     plot(deg0int,gamma(4,:,1),'b-', 'linewidth', 2);
%     legend('VV','HV','VH','HH');
%     xlabel('\theta_{obs} (\circ)');
%     ylabel('\gamma (dB)');
%     title('biscatic scattering on incidence plane');
%     set(gca,'XMinorTick','on');
%     set(gca,'YMinorTick','on');
%     xlim([-80 80]); grid on; box on;
% 
% $Data: July 4, 2013
% $Revision: 1.0
% 
    

%% debug - preset parameters
% ============== debgu begin
% deg0inc = 40;
% freq = 9.6;
% depth = [1 10 100]';
% dia = 0.3;
% fractv = 0.2;
% perms = 3.15 + 0.001i;
% tau = 0.15;
% rsfile = 'roughsurface_X.mat';
% ndeg = 16;
% Mmax = 4;
% debug end =================

%% interprete the input variable
bFlat = isempty(rsfile);
if (~bFlat) && isempty(perm_gnd)
    perm_gnd = 3.2 + 1i*0.002;     % do not take effect with rsbc
end

[nd,nl] = size(depth);
[nl_c, ~] = size(P11_a);
if nl ~= length(perm_eff_a) || nl ~= nl_c
    error('Layer number mis-match');
end

% for off diagonal elements of phase matrix.
if nargin < 17
    np = length(bmua);
    P13_a = zeros(nl,np); 
    P14_a = zeros(nl,np); 
    P23_a = zeros(nl,np); 
    P24_a = zeros(nl,np); 
    P31_a = zeros(nl,np); 
    P32_a = zeros(nl,np); 
    P41_a = zeros(nl,np); 
    P42_a = zeros(nl,np); 
end

% dia = App(dia,nl);
% fractv = App(fractv,nl);
% perms = App(perms,nl);
% tau = App(tau,nl);

% locate output
sigma = zeros(4,nd);
gamma = zeros(4,ndeg*2,nd);
%% depth independent variable

%% phase matrix preparation
% dependent parameters:
%   nl - number of layers
%   freq - frequency
%   dia, fractv, perms, tau - grain parameters, all arrays

% fprintf('Calculate Phase Matrix ...\n');
% Nmu = 361; % number of angles in calculating phase matrix.
% P11_a = zeros(nl,Nmu);
% P22_a = zeros(nl,Nmu);
% P33_a = zeros(nl,Nmu);
% P44_a = zeros(nl,Nmu);
% P34_a = zeros(nl,Nmu);
% P43_a = zeros(nl,Nmu);
% P12_a = zeros(nl,Nmu); % - added
% P21_a = zeros(nl,Nmu); % - added
% ka = zeros(nl,1);
% ke = zeros(nl,1);
% ks = zeros(nl,1);
% keff_a = zeros(nl,1);
% 
% for in = 1:nl
%     [keff_a(in),ka(in),ks(in),ke(in),...
%         bmua,P11_a(in,:),P22_a(in,:),...
%         P33_a(in,:),P44_a(in,:),P34_a(in,:),P43_a(in,:)] = ...
%         qcaphase_sticky_bigP12_new(freq,dia(in),fractv(in),perms(in),tau(in),Nmu);
% end

%% extract effective medium parameters
% lambda      = 30/freq;          % wavelength in centimeter;
% k          =2*pi/lambda;     %wave number
% Kprime     =real(keff_a); 
% % epsilon_drysnow = (Kprime/k).^2;                  %permittivity of dry snow
% % eps1 = epsilon_drysnow;
% % epsilon_2 = 3.2 + 1i*0.002;     % do not take effect with rsbc
% k2 = sqrt(epsilon_2)*k;
% ke = ka + ks;
% % albedo = ks ./ ke;
% Kp = Kprime ;    % Kp =Kprime1,Kprime2,...Kprimen
% % K = Kp+ 1i * ke / 2;
% kp = [ k; Kp; k2];  % kp=k0,k1,k2,...kn,kg

% use normalized frequency domain quantitiy.
k = 1;
Kp = sqrt(real(perm_eff_a))*k;
k2 = sqrt(perm_gnd)*k;
kp = [k;Kp;k2];
ke = ka + ks;
%% snell's law to determine local incident angles
% new dependent paras:
%   deg0inc - incident angle

tai0inc = deg0inc * pi / 180;
taiinc = asin ( k * sin ( tai0inc ) ./ Kp );
ctinc = cos(taiinc);
% ct0inc = cos(tai0inc);

muinc=zeros(2,nl);
muinc(1,:)=cos(taiinc);
muinc(2,:)=-muinc(1,:);
%% Fourier Transform of Phase Matrx, Discrete Coordinate Quadrature
fprintf('Fourier Series of Phase Matrix ...\n');
% new dependent paras:
%   ndeg - half of the number of quadrature angles
%   Mmax - Fourier Series Order

nquad=ndeg*2;
ndeg2=ndeg*2;

% [x,w]=GAULEG(-1,1,nquad);
[x,w] = GLNodeWt(nquad);
% renumber as follows
%first half is upward going  and from small angle to large angle
%second half downward but in order of angle as the first half
mut=x(ndeg2:-1:(ndeg+1));
aat=w(ndeg2:-1:(ndeg+1));
mub=x(1:ndeg);
aab=w(1:ndeg);
mu=[mut;mub];
aa=[aat;aab];

P0_a = zeros(ndeg2,ndeg2,4,4,nl);
PC_a = zeros(Mmax,ndeg2,ndeg2,4,4,nl);
PS_a = zeros(Mmax,ndeg2,ndeg2,4,4,nl);

ndeginc = 2; % taiinc, and pi - taiinc
P0inc_a = zeros(ndeg2,ndeginc,4,4,nl);
PCinc_a = zeros(Mmax,ndeg2,ndeginc,4,4,nl);
PSinc_a = zeros(Mmax,ndeg2,ndeginc,4,4,nl);

%number of phi angle when transform from 12 frame to principle frame.
% nphs = 61; 
nphs = 181; 
for jn=1:nl
    [P0_a(:,:,:,:,jn),PC_a(:,:,:,:,:,jn),PS_a(:,:,:,:,:,jn)] = ...
        transf12p_fc_all(bmua,P11_a(jn,:),P22_a(jn,:),P12_a(jn,:),P21_a(jn,:),P33_a(jn,:),P34_a(jn,:),mu,mu,nphs,Mmax,...
        P13_a(jn,:), P14_a(jn,:), P23_a(jn,:), P24_a(jn,:), P31_a(jn,:), P32_a(jn,:), P41_a(jn,:), P42_a(jn,:));
    P0_a(:,:,:,:,jn) = P0_a(:,:,:,:,jn)*ks(jn);
    PC_a(:,:,:,:,:,jn) = PC_a(:,:,:,:,:,jn)*ks(jn);
    PS_a(:,:,:,:,:,jn) = PS_a(:,:,:,:,:,jn)*ks(jn);
    
    [P0inc_a(:,:,:,:,jn),PCinc_a(:,:,:,:,:,jn),PSinc_a(:,:,:,:,:,jn)]...
        = transf12p_fc_all(bmua,P11_a(jn,:),P22_a(jn,:),P12_a(jn,:),P21_a(jn,:),P33_a(jn,:),P34_a(jn,:),[ctinc(jn) -ctinc(jn)],mu,nphs,Mmax,...
        P13_a(jn,:), P14_a(jn,:), P23_a(jn,:), P24_a(jn,:), P31_a(jn,:), P32_a(jn,:), P41_a(jn,:), P42_a(jn,:));  
    P0inc_a(:,:,:,:,jn) = P0inc_a(:,:,:,:,jn)*ks(jn);
    PCinc_a(:,:,:,:,:,jn) = PCinc_a(:,:,:,:,:,jn)*ks(jn);
    PSinc_a(:,:,:,:,:,jn) = PSinc_a(:,:,:,:,:,jn)*ks(jn);
end

%% prepare reflection coefficient on quadrature angles, interpolation
% on flat boundaries
[rup,rdown,tup,tdown,rup0,rdown0,tup0,tdown0] = ref_tran(kp,mut,nl,ndeg);
% rup dimensions: (8*ndeg,8*ndeg,nl), sequence: r10,r21,...rn,n-1,

% pre-calculate interpolation matrix
% [miu,~] = GLNodeWt(2*ndeg);
[Sd,Su,Sd0,Su0] = spline_stokes(nl,mut,Kp,k);
% Sd dimensions: (8*ndeg,8*ndeg,nl-1), sequence of Sd: S12,S23,...Sn-1,n Su:S21,S32,...Sn,n-1

%% modification due to rough surface, reflection, and inc-biscatic scat
if ~bFlat
    load(rsfile);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % modify rdown0(:,:,nl) and rdown(:,:,nl) to coherent rough surface
    % reflection coefficient
    tagd = acos(mut); % quadrature angles (0, 90)
    [rdown0(:,:,nl),rdown(:,:,nl)] = ref_coh(rc_rs,tai_rs,tagd);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prepare coherent scattering coefficient for the rough snow-ground
    % surface for the relating to the incident angle
    r12inc = ref_coh_inc(rc_rs,tai_rs,taiinc(nl));

    % prepare incoherent bistatic scattering coefficient for the rough
    % surface snow-ground interface, Fourier series expansion
    fprintf('Fourier Transform of Rough Surface Bistatic Scattering ...\n');
    tagd2 = acos((mu + 1)/2); % quadrature angles (0, 90)
    [gammaic0t, gammaic0tinc, gammaict, gammaictinc] = ...
    transf_incoh(gammainc_rs,tas_rs,phs_rs,tai_rs,tagd2,taiinc(end),Mmax);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Calculate Matrix Elements of Eigenvalue problem
fprintf('Calculate Matrix Elements of Eigenvalue Problem ...\n');
Am_a = zeros(16*ndeg,16*ndeg,Mmax,nl);
Am0_a = zeros(8*ndeg,8*ndeg,nl);
for jn = 1:nl
    [Am0_a(:,:,jn), Am_a(:,:,:,jn)] = ...
        CalcAm(mu,aa,P0_a(:,:,:,:,jn),PC_a(:,:,:,:,:,jn),PS_a(:,:,:,:,:,jn));
end

%% preliminary job for the eigenvalue problem, particular problem, and
% boundary condition problem
fprintf('Solving Eigenvalue Problem ...\n');
% eigenvalue and eigen vector
eigvp0=zeros(4*ndeg,nl);
eigvm0=zeros(4*ndeg,nl);
vecup0=zeros(4*ndeg,4*ndeg,nl);
vecdp0=zeros(4*ndeg,4*ndeg,nl);
vecum0=zeros(4*ndeg,4*ndeg,nl);
vecdm0=zeros(4*ndeg,4*ndeg,nl);

eigvp=zeros(4*ndeg2,nl,Mmax);
eigvm=zeros(4*ndeg2,nl,Mmax);
vecup=zeros(4*ndeg2,4*ndeg2,nl,Mmax);
vecdp=zeros(4*ndeg2,4*ndeg2,nl,Mmax);
vecum=zeros(4*ndeg2,4*ndeg2,nl,Mmax);
vecdm=zeros(4*ndeg2,4*ndeg2,nl,Mmax);

% particular problem
AP10_a = zeros(4*ndeg2,4*ndeg2,nl); % M - source term
AP20_a = zeros(4*ndeg2,4*ndeg2,nl); % N - source term
AP1_a = zeros(8*ndeg2,8*ndeg2,nl,Mmax); % M - source term
AP2_a = zeros(8*ndeg2,8*ndeg2,nl,Mmax); % N - source term

for jn = 1:nl
    % zeroth harmonic
    [eigvp0(:,jn), vecup0(:,:,jn), vecdp0(:,:,jn),...
            eigvm0(:,jn), vecum0(:,:,jn), vecdm0(:,:,jn),...
            AP10_a(:,:,jn), AP20_a(:,:,jn)] ...
            = EP(Am0_a(:,:,jn),ke(jn),mu,muinc(:,jn));
    % m-th harmonic
    for m = 1:Mmax
        [eigvp(:,jn,m), vecup(:,:,jn,m), vecdp(:,:,jn,m),...
            eigvm(:,jn,m), vecum(:,:,jn,m), vecdm(:,:,jn,m),...
            AP1_a(:,:,jn,m), AP2_a(:,:,jn,m)] ...
            = EP(Am_a(:,:,m,jn),ke(jn),mu,muinc(:,jn));
    end
end


% boundary condition problem
% top and bottom boundaries
fprintf('Prepare for Boundary Problem ...\n');
% A10_a = zeros(4*ndeg,8*ndeg);
% An0_a = zeros(4*ndeg,8*ndeg);
A1_a = zeros(8*ndeg,16*ndeg,Mmax);
An_a = zeros(8*ndeg,16*ndeg,Mmax);

% zero-th harmonic
[A10_a, An0_a] = BPext(vecup0(:,:,1),vecdp0(:,:,1),vecum0(:,:,1),vecdm0(:,:,1),rup0(:,:,1),...
    vecup0(:,:,nl),vecdp0(:,:,nl),vecum0(:,:,nl),vecdm0(:,:,nl),rdown0(:,:,nl));
if ~bFlat
    % Modify An due to rough snow-ground boundary.
    As0 = BPrs(gammaic0t,mu,aa,vecdm0(:,:,nl),vecdp0(:,:,nl));
    An0_a = An0_a - As0;
end

% m-th harmonic
for m = 1:Mmax
    [A1_a(:,:,m), An_a(:,:,m)] = BPext(vecup(:,:,1,m),vecdp(:,:,1,m),vecum(:,:,1,m),vecdm(:,:,1,m),rup(:,:,1),...
        vecup(:,:,nl,m),vecdp(:,:,nl,m),vecum(:,:,nl,m),vecdm(:,:,nl,m),rdown(:,:,nl));
    if ~bFlat
        % Modify An due to rough snow-ground boundary.
        As = BPrs(gammaict(:,:,m),mu,aa,vecdm(:,:,nl,m),vecdp(:,:,nl,m));
        An_a(:,:,m) = An_a(:,:,m) - As;
    end
end

% internal boundaries
Au0_a=zeros(8*ndeg,8*ndeg,nl-1);
Ad0_a=zeros(8*ndeg,8*ndeg,nl-1);
Au_a=zeros(16*ndeg,16*ndeg,nl-1,Mmax);
Ad_a=zeros(16*ndeg,16*ndeg,nl-1,Mmax);

for jn = 1:nl - 1
    % zeroth harmonic
    [Au0_a(:,:,jn), Ad0_a(:,:,jn)] = BPint(vecup0(:,:,jn),vecdp0(:,:,jn),...
        vecum0(:,:,jn),vecdm0(:,:,jn),rdown0(:,:,jn),tdown0(:,:,jn),Sd0(:,:,jn),...
        vecup0(:,:,jn + 1),vecdp0(:,:,jn + 1),vecum0(:,:,jn + 1),vecdm0(:,:,jn + 1),...
        rup0(:,:,jn + 1),tup0(:,:,jn + 1),Su0(:,:,jn));
    % m-th harmonic
    for m = 1:Mmax
        [Au_a(:,:,jn,m), Ad_a(:,:,jn,m)] = BPint(vecup(:,:,jn,m),vecdp(:,:,jn,m),...
            vecum(:,:,jn,m),vecdm(:,:,jn,m),rdown(:,:,jn),tdown(:,:,jn),Sd(:,:,jn),...
            vecup(:,:,jn + 1,m),vecdp(:,:,jn + 1,m),vecum(:,:,jn + 1,m),vecdm(:,:,jn + 1,m),...
            rup(:,:,jn + 1),tup(:,:,jn + 1),Su(:,:,jn));
    end
end


%% following are depth dependent procedures.
for id = 1:nd
d = depth(id,:);
fprintf('depth = %.1f\n',sum(d));
%% zero-order, calculate M, N
% new dependent parameters:
%   d - snow depth, vector configuring the layer.
fprintf('\tSolving reduced intensity ...\n');
if bFlat
    [a,b,c,Hj] = abc(nl,ke,kp,taiinc,d);
else
    [a,b,c,Hj] = abc_rs(nl,ke,kp,taiinc,d,r12inc);
end
A0 = ( b * c - a ) \ b;
% A0 = b / ( b * c - a );
B0 = - b \ a * A0;
% M_i, N_i are of dimensions (4,nl,2)
[M_i,N_i] = MN(nl,A0,B0,k,Kp,taiinc,tai0inc,Hj);

%% calculate RHS elements for particular problem
% depending on M and N vector, and thus on depth
RHS1_a = zeros(16*ndeg,2,Mmax,nl);
RHS2_a = zeros(16*ndeg,2,Mmax,nl);
RHS10_a=zeros(8*ndeg,2,nl);
RHS20_a=zeros(8*ndeg,2,nl);

for jn = 1:nl
    M0inc = reshape(M_i(:,jn,:),4,2);
    N0inc = reshape(N_i(:,jn,:),4,2);
    [rhs1,rhs2,rhs10,rhs20] = CalcRHS(mu,P0inc_a(:,:,:,:,jn),...
        PCinc_a(:,:,:,:,:,jn),PSinc_a(:,:,:,:,:,jn),M0inc,N0inc);
    RHS1_a(:,:,:,jn) = rhs1;
    RHS2_a(:,:,:,jn) = rhs2;
    RHS10_a(:,:,jn) = rhs10;
    RHS20_a(:,:,jn) = rhs20;
end

%% Seek particular solution, and solve boundary condition problem

%% zero-th harmonic
fprintf('\tZero-th harmonic Particular Sol. and BVP...\n');
dimm = 4*ndeg2;
dimh = dimm/2;
rS = 1:dimh;
rB = dimh + 1:dimm;

% upward specific intensity at the top boundary in the air,
% first dimmension, number of angles,
% second dimension, for v and h pol,
% I1U0 = zeros(dimh,2);

% seek particular solution.
part1u0 = zeros(dimh,2,nl);  % the second dimension is for v and h pol.
part1d0 = zeros(dimh,2,nl);
part2u0 = zeros(dimh,2,nl);
part2d0 = zeros(dimh,2,nl);

for nk = 1:nl
    % the upward downward sort is from mu factor in RHS.
    part10 = AP10_a(:,:,nk)\RHS10_a(:,:,nk);
    part1u0(:,:,nk) = part10(rS,:);
    part1d0(:,:,nk) = part10(rB,:);

    part20 = AP20_a(:,:,nk)\RHS20_a(:,:,nk);
    part2u0(:,:,nk) = part20(rS,:);
    part2d0(:,:,nk) = part20(rB,:);
end

% solving simulaneous equation by imposing boundary conditions
% Fill the BC element
[Au0,Ad0,B0,A10,An0,B10,Bn0] = FillBC(d,eigvp0,eigvm0,...
    Au0_a,Ad0_a,A10_a,An0_a,...
    ke,ctinc,part1d0,part1u0,part2d0,part2u0,...
    rdown0,rup0,tdown0,tup0,Sd0,Su0);
if ~bFlat
    % modify Bn due to rough surface
    Bs0 = BCrs(ke(nl),d(nl),ctinc(nl),...
        part1d0(:,:,nl),part2d0(:,:,nl),reshape(N_i(:,nl,:),4,2),...
        mu,aa,gammaic0t,gammaic0tinc);
    Bn0 = Bn0 + Bs0;
end
% solve the boundary problem
coea0 = SolveBC(Au0,Ad0,B0,A10,An0,B10,Bn0);
% construct the specific intensity at the top layer or in the air
I1U0 = CalcIU(coea0(1:dimm,:),vecup0(:,:,1),vecum0(:,:,1),eigvm0(:,1),...
    ke(1),d(1),ctinc(1),part1u0(:,:,1),part2u0(:,:,1),tup0(:,:,1));


%% m-th harmonic
dimm=8*ndeg2;
dimh=dimm/2;

rS = 1:dimh;
rB = dimh + 1:dimm;

% upward specific intensity at the top boundary in the air,
% first dimmension, number of angles,
% second dimension, for v and h pol,
% third dimension, for the order of harmonics.
I1U = zeros(dimh,2,Mmax); 

for m = 1:Mmax
    fprintf('\t%d-th harmonic Particular Sol. and BVP...\n',m);
    % seek particular solution.
    part1u = zeros(dimh,2,nl);  % the second dimension is for v and h pol.
    part1d = zeros(dimh,2,nl);
    part2u = zeros(dimh,2,nl);
    part2d = zeros(dimh,2,nl);
    
    for nk = 1:nl
        % the upward downward sort is from mu factor in RHS.
        part1 = AP1_a(:,:,nk,m)\RHS1_a(:,:,m,nk);
        part1u(:,:,nk) = part1(rS,:);
        part1d(:,:,nk) = part1(rB,:);
        
        part2 = AP2_a(:,:,nk,m)\RHS2_a(:,:,m,nk);
        part2u(:,:,nk) = part2(rS,:);
        part2d(:,:,nk) = part2(rB,:);
    end
    
    % solving simulaneous equation by imposing boundary conditions
    % Fill the BC element
    [Au,Ad,B,A1,An,B1,Bn] = FillBC(d,eigvp(:,:,m),eigvm(:,:,m),...
        Au_a(:,:,:,m),Ad_a(:,:,:,m),A1_a(:,:,m),An_a(:,:,m),...
        ke,ctinc,part1d,part1u,part2d,part2u,...
        rdown,rup,tdown,tup,Sd,Su);
    if ~bFlat
        % modify Bn due to rough surface
        Bs = BCrs(ke(nl),d(nl),ctinc(nl),...
            part1d(:,:,nl),part2d(:,:,nl),reshape(N_i(:,nl,:),4,2),...
            mu,aa,gammaict(:,:,m),gammaictinc(:,:,m));
        Bn = Bn + Bs;
    end
    % solve the boundary problem
    coeam = SolveBC(Au,Ad,B,A1,An,B1,Bn);
    % construct the specific intensity at the top layer or in the air
    I1U(:,:,m) = CalcIU(coeam(1:dimm,:),vecup(:,:,1,m),vecum(:,:,1,m),eigvm(:,1,m),...
        ke(1),d(1),ctinc(1),part1u(:,:,1),part2u(:,:,1),tup(:,:,1));
    
end

%% organize the results
fprintf('\tTransform to incident plane...\n');
% Transfrom the specific intensity spectrum into intensity in the 
% incident plane (phi = pi, theta[0, pi/2]) (phi = 0, theta [0,pi/2])
% by summing up Fourier Series
% Relate the angle in region 1 and zero by Snell's law, considering
% critical angle effect.
% Also calculate backscattering by interpolation.
[sigma_,gamma_,deg0int] = TransIncPlane(I1U0,I1U,mut,Kp(1)/k,tai0inc);
% record the results
sigma(1,id) = sigma_(1,1);
sigma(2,id) = sigma_(2,1);
sigma(3,id) = sigma_(1,2);
sigma(4,id) = sigma_(2,2);
ndeg0int = length(deg0int);
gamma(1,1:ndeg0int,id) = reshape(gamma_(1,1,:),1,ndeg0int,1);
gamma(2,1:ndeg0int,id) = reshape(gamma_(2,1,:),1,ndeg0int,1);
gamma(3,1:ndeg0int,id) = reshape(gamma_(1,2,:),1,ndeg0int,1);
gamma(4,1:ndeg0int,id) = reshape(gamma_(2,2,:),1,ndeg0int,1);
end

gamma = gamma(:,1:ndeg0int,:);
fprintf('Finished.\n');
end