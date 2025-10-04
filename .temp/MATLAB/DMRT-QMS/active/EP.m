function [eigvp, vecup, vecdp, eigvm, vecum, vecdm, AP1, AP2] ...
    = EP(Am,ke,mu,muinc)
% preliminary job for the eigenvalue problem, particular problem

% eigenvalue problem ===============================================
[dimm, ~] = size(Am); % 8*ndeg2
ndeg2 = length(mu);
% dimm = 8*ndeg2;
% dimh = 4*ndeg2;
dimh = dimm/2;
Nel = dimm/ndeg2;

% [dimm,~] = size(Am); % 8*ndeg2
% dimh = dimm/2;       % 4*ndeg2
% ndeg2 = dimm/8;     % number of discrete angles
% nc = 2; % 1-cosine; 2-sine;

eigvp=zeros(dimh,1);
eigvm=zeros(dimh,1);
vecup=zeros(dimh,dimh);
vecdp=zeros(dimh,dimh);
vecum=zeros(dimh,dimh);
vecdm=zeros(dimh,dimh);

% EA is the matrix for the eigenvalue problem
%, it is A with chnages in the diagonal elements
sectheta = folder(1./mu,Nel);
secthe=diag(sectheta);
EA = Am - ke*secthe;
[vec value]=eig(EA);

% sort eigenvalue and vector
%8 N postive eigenvalues, for each positive eigenvalue,we have
%a single 16N x1 eigenvector, the top half is upward and the
%bottom half is down ward Thus vecup is 8N x 8N 8N negative
%eigenvalues eigenvector of upward going with positive
%eigenvalue,  vecup eigenvector downward,  positive eigenvalue
%vecdp eigenvector of upward negative (minus) eigenvalue vecum
%eiegenvector downward,  negative eigenvalue vecdm
icolp=0;
icolm=0;
for icol=1:dimm  %2345
    q=value(icol,icol);
    rq=real(q);
    aq=abs(rq); 
    iq=imag(q);
    if  ((aq>=1e-8)&&(rq>0))||((aq<1e-8)&&(iq>0)) %112
        icolp=icolp+1;
        eigvp(icolp)=q;
        vecup(:,icolp)=vec(1:dimh,icol);
        vecdp(:,icolp)=vec(1+dimh:2*dimh,icol);
    else
        icolm=icolm+1;
        eigvm(icolm)=q;
        vecum(:,icolm)=vec(1:dimh,icol);
        vecdm(:,icolm)=vec(1+dimh:2*dimh,icol);
    end %112
end %2345

% particular problem ===========================================
% particular solutions find the particular solutions.  The
% matrix is the same as EA with diagonal elements modified
AP1 = - EA - ke/muinc(1)*eye(dimm);  % M - source term
AP2 = - EA - ke/muinc(2)*eye(dimm);  % N - source term

end

function [y] = folder(x,s)
y = reshape(repmat((x(:))',s,1),length(x)*s,1);
end
