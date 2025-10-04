function [k0,keff_a,ka,ks,bmua,P11_a,P22_a,P33_a,P34_a,P12_a,P21_a] ...
    = CalcPM(fGHz, dia, fractv, perms, tau, Nmu)
% Calculate phase matrix with sticky qca approach.
% The phase matrix are normalized by kappas.
% keff_a, ka and ks are in (1/cm)

fprintf('Calculate Phase Matrix ...\n');
if nargin < 6
    Nmu = 65; % number of angles in calculating phase matrix.
end

nl = max([length(dia),length(fractv),length(perms),length(tau)]);
dia = App(dia,nl);
fractv = App(fractv,nl);
perms = App(perms,nl);
tau = App(tau,nl);

P11_a = zeros(nl,Nmu);
P22_a = zeros(nl,Nmu);
P33_a = zeros(nl,Nmu);
% P44_a = zeros(nl,Nmu); % P44 = P33;
P34_a = zeros(nl,Nmu);
% P43_a = zeros(nl,Nmu); % P43 = -P34;
P12_a = zeros(nl,Nmu); % - added
P21_a = zeros(nl,Nmu); % - added
ka = zeros(nl,1);
ke = zeros(nl,1);
ks = zeros(nl,1);
keff_a = zeros(nl,1);

lambda = 30/fGHz; % cm
k0 = 2*pi/lambda; % 1/cm
for in = 1:nl
    [keff_a(in),ka(in),ks(in),ke(in),...
        bmua,P11,P22,P33,~,P34,~] = ...
        qcaphase_sticky_bigP12_new(fGHz,dia(in),fractv(in),perms(in),tau(in),Nmu);
    
    P11_a(in,:) = P11/ks(in);
    P22_a(in,:) = P22/ks(in);
    P33_a(in,:) = P33/ks(in);
    P34_a(in,:) = P34/ks(in);
end


end