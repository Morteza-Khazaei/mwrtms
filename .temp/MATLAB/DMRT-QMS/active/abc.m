function [a,b,c,Hj]=abc(nl,ke,kp,taiinc,d)

h=zeros(8,8,nl-1); %h21, h32, ...,hn,n-1
k = kp(1);
Kp = kp(2:nl+1);
% every m
%Sd:S12,S23,...Sn-1,n Su:S21,S32,...Sn,n-1
H=eye(8); % hn,n-1*hn-1,n-2...*h21
Hj=zeros(8,8,nl-1); %h21,h32*h21,h43*h32*h21,...hnn-1*...*h21
for in=1:nl-1
    h(1:4,1:4,in)=exp(ke(in+1)*sec(taiinc(in+1))*d(in+1))*inv(tranmatrix(Kp(in+1),Kp(in),taiinc(in+1)));
    h(1:4,5:8,in)=-exp(ke(in+1)*sec(taiinc(in+1))*d(in+1)-ke(in)*sec(taiinc(in))*d(in))*refmatrix(Kp(in),Kp(in+1),taiinc(in))*inv(tranmatrix(Kp(in+1),Kp(in),taiinc(in+1)));
    h(5:8,1:4,in)=refmatrix(Kp(in+1),Kp(in),taiinc(in+1))*inv(tranmatrix(Kp(in+1),Kp(in),taiinc(in+1)));
    h(5:8,5:8,in)=(tranmatrix(Kp(in),Kp(in+1),taiinc(in))-refmatrix(Kp(in+1),Kp(in),taiinc(in+1))*refmatrix(Kp(in),Kp(in+1),taiinc(in))*inv(tranmatrix(Kp(in+1),Kp(in),taiinc(in+1))))*exp(-ke(in)*sec(taiinc(in))*d(in));
    H=h(:,:,in)*H;
    Hj(:,:,in)=H;
end
H11=H(1:4,1:4);
H12=H(1:4,5:8);
H21=H(5:8,1:4);
H22=H(5:8,5:8);
 % kp=k0,k1,k2,...kn,kg
a=H11-exp(-ke(nl)*sec(taiinc(nl))*d(nl))*refmatrix(kp(nl+1),kp(nl+2),taiinc(nl))*H21;
b=H12-exp(-ke(nl)*sec(taiinc(nl))*d(nl))*refmatrix(kp(nl+1),kp(nl+2),taiinc(nl))*H22;
c=-exp(-ke(1)*sec(taiinc(1))*d(1))*refmatrix(Kp(1),k,taiinc(1));