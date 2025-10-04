function [M_i,N_i]=MN(nl,A0,B0,k,Kp,taiinc,tai0inc,Hj)

Mi=zeros(4,4,nl);
Ni=zeros(4,4,nl);
M_i=zeros(4,nl,2);
N_i=zeros(4,nl,2);


%I0 is incident Stokes vector
for ipol=1:2  %1000
    if  (ipol==1)
        I0=[1;0;0;0];
    end
    if  (ipol==2)
        I0=[0;1;0;0];
    end
    Mi(:,:,1) = A0;
    M_i(:,1,ipol) = A0 * tranmatrix( k ,Kp(1),tai0inc ) * I0 * ( k / Kp(1) )^2 * cos(tai0inc) / cos( taiinc(1) );
    Ni(:,:,1) = B0;
    N_i(:,1,ipol) = B0 * tranmatrix( k , Kp(1),tai0inc ) * I0 * ( k / Kp(1) )^2 * cos(tai0inc) / cos( taiinc(1) );
    for in=1:nl-1
        Mi(:,:,in+1) = Hj(1:4,1:4,in) * A0 + Hj(1:4,5:8,in) * B0;
        M_i(:,in+1,ipol) = Mi(:,:,in+1) * tranmatrix( k , Kp(1) , tai0inc ) * I0 * ( k / Kp (in+1) )^2 * cos(tai0inc) / cos( taiinc(in+1) );
        Ni(:,:,in+1) = Hj(5:8,1:4,in) * A0 + Hj(5:8,5:8,in) * B0;
        N_i(:,in+1,ipol) = Ni(:,:,in+1) * tranmatrix( k , Kp(1) ,tai0inc ) * I0 * (k/Kp(in+1))^2 * cos(tai0inc) / cos( taiinc(in+1) );
    end
end   %1000