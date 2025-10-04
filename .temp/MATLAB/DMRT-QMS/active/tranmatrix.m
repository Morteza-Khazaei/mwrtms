%tranmatrix.m
% transmissivity matrix for wave incident at angle tai from kinc medium to
% ktran medium.  tai is measured in kinc medium

function [t]=tranmatrix(kinc,ktran,tai)
%kinc is real and ktran is real for transmission to be meaninful and Snell's law
% this is unlike reflection.  In transmission, the transmitted medium must have real permitivity
% or close to real permittivity

         t=zeros(4);
   kiz=kinc*cos(tai);
    krho=kinc*sin(tai);
    thetac=pi/2;
    %find critical angle
    if (kinc>ktran)  %101
        thetac=asin(ktran/kinc);
    end  %101
    if  (tai<thetac)  %300
        tait=asin(kinc*sin(tai)/ktran);
    ktranz=sqrt(ktran^2-krho*krho);
        if (imag(ktranz)<0) %101
        ktranz=-ktranz;
         end %101
    RV=(ktran*ktran*kiz-kinc*kinc*ktranz)/(ktran*ktran*kiz+kinc*kinc*ktranz);
    RH=(kiz-ktranz)/(kiz+ktranz);
     t(1,1)=1-abs(RV)^2;
     t(2,2)=1-abs(RH)^2;
     TV=1+RV;
     TH=1+RH;
     cfac=cos(tait)/cos(tai);
     t(3,3)=cfac*real(TV*conj(TH));
     t(4,4)=t(3,3);
     t(3,4)=-cfac*imag(TV*conj(TH));
     t(4,3)=-t(3,4);
     %divergencde of beam factor
     t=t*ktran^2/kinc^2;
 end  % 300