function [RHS1,RHS2,RHS10,RHS20] = CalcRHS(mu,P0inc,PCinc,PSinc,M0inc,N0inc)

[Mmax,ntas,~,~,~] = size(PCinc);
ndeg = ntas/2;

RHS1=zeros(16*ndeg,2,Mmax);
% RHS1v=zeros(16*ndeg,Mmax);
% RHS1h=zeros(16*ndeg,Mmax);
RHS2=zeros(16*ndeg,2,Mmax);
% RHS2v=zeros(16*ndeg,Mmax);
% RHS2h=zeros(16*ndeg,Mmax);

RHS10=zeros(8*ndeg,2);
RHS20=zeros(8*ndeg,2);
% RHS1v0=zeros(8*ndeg,1);
% RHS1h0=zeros(8*ndeg,1);
% RHS2v0=zeros(8*ndeg,1);
% RHS2h0=zeros(8*ndeg,1);

for idegs=1:ntas  %201
    ids_ = block(idegs,4);
    RHS10(ids_,:)=reshape(P0inc(idegs,1,:,:),4,4)*M0inc/mu(idegs);
    RHS20(ids_,:)=reshape(P0inc(idegs,2,:,:),4,4)*N0inc/mu(idegs);
    
%     RHS1v0(ids_)=P0inc(idegs,1,:,:)*M0incv/mu(idegs);
%     RHS1h0(ids_)=P0inc(idegs,1,:,:)*M0inch/mu(idegs);
%     
%     RHS2v0(ids_)=P0inc(idegs,2,:,:)*N0incv/mu(idegs);
%     RHS2h0(ids_)=P0inc(idegs,2,:,:)*N0inch/mu(idegs);
    
    ids = block(idegs,8);
    for m = 1:Mmax
        RHS1(ids(1:4),:, m)=reshape(PCinc(m,idegs,1,:,:),4,4)*M0inc/mu(idegs);
        RHS1(ids(5:8),:, m)=reshape(PSinc(m,idegs,1,:,:),4,4)*M0inc/mu(idegs);
        RHS2(ids(1:4),:, m)=reshape(PCinc(m,idegs,2,:,:),4,4)*N0inc/mu(idegs);
        RHS2(ids(5:8),:, m)=reshape(PSinc(m,idegs,2,:,:),4,4)*N0inc/mu(idegs);
        
%         RHS1v(ids(1:4), m)=PCinc(m,idegs,1,:,:)*M0incv/mu(idegs);
%         RHS1v(ids(5:8), m)=PSinc(m,idegs,1,:,:)*M0incv/mu(idegs);
%         RHS1h(ids(1:4), m)=PCinc(m,idegs,1,:,:)*M0inch/mu(idegs);
%         RHS1h(ids(5:8), m)=PSinc(m,idegs,1,:,:)*M0inch/mu(idegs);
%         
%         RHS2v(ids(1:4), m)=PCinc(m,idegs,2,:,:)*N0incv/mu(idegs);
%         RHS2v(ids(5:8), m)=PSinc(m,idegs,2,:,:)*N0incv/mu(idegs);
%         RHS2h(ids(1:4), m)=PCinc(m,idegs,2,:,:)*N0inch/mu(idegs);
%         RHS2h(ids(5:8), m)=PSinc(m,idegs,2,:,:)*N0inch/mu(idegs);
    end
end
end

function id = block(x,s)
id = (x-1)*s + 1:x*s;
end