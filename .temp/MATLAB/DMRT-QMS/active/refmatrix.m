%refmatrix.m
%reflectivity matrx
%kinc is real and ktran can be complex
% reflectivity matrix for wave incident at angle tai from kinc medium to
% ktran medium.  tai is measured in kinc medium

function [r]=refmatrix(kinc,ktran,tai)
    kiz=kinc*cos(tai);
    krho=kinc*sin(tai);
    ktranz=sqrt(ktran^2-krho*krho);
    if (imag(ktranz)<0)
        ktranz=-ktranz;
    end
    RV=(ktran*ktran*kiz-kinc*kinc*ktranz)/(ktran*ktran*kiz+kinc*kinc*ktranz);
    RH=(kiz-ktranz)/(kiz+ktranz);
     r=zeros(4);
     r(1,1)=abs(RV)^2;
     r(2,2)=abs(RH)^2;
     r(3,3)=real(RV*conj(RH));
     r(4,4)=r(3,3);
     r(3,4)=-imag(RV*conj(RH));
     r(4,3)=-r(3,4);