function [rup,rdown,tup,tdown,rup0,rdown0,tup0,tdown0]=ref_tran(kp,mu,nl,ndeg)

rup = zeros(8*ndeg,8*ndeg,nl);  %r10,r21,...rn,n-1,
rdown = zeros(8*ndeg,8*ndeg,nl);  %r12,r23,...rn,n+1
tup = zeros(8*ndeg,8*ndeg,nl);  %t10,t21,...tn,n-1,
tdown = zeros(8*ndeg,8*ndeg,nl);  %t12,t23,...tn,n+1
% kp=k0,k1,k2,...kn
for in=1:nl
    for ideg=1:ndeg
        tai=acos(mu(ideg));
        rup(1+8*(ideg-1):4+8*(ideg-1),1+8*(ideg-1):4+8*(ideg-1),in)=refmatrix(kp(in+1),kp(in),tai);
        rup(5+8*(ideg-1):8+8*(ideg-1),5+8*(ideg-1):8+8*(ideg-1),in)=refmatrix(kp(in+1),kp(in),tai);
        rdown(1+8*(ideg-1):4+8*(ideg-1),1+8*(ideg-1):4+8*(ideg-1),in)=refmatrix(kp(in+1),kp(in+2),tai);
        rdown(5+8*(ideg-1):8+8*(ideg-1),5+8*(ideg-1):8+8*(ideg-1),in)=refmatrix(kp(in+1),kp(in+2),tai);
        
        tup(1+8*(ideg-1):4+8*(ideg-1),1+8*(ideg-1):4+8*(ideg-1),in)=tranmatrix(kp(in+1),kp(in),tai);
        tup(5+8*(ideg-1):8+8*(ideg-1),5+8*(ideg-1):8+8*(ideg-1),in)=tranmatrix(kp(in+1),kp(in),tai);
        tdown(1+8*(ideg-1):4+8*(ideg-1),1+8*(ideg-1):4+8*(ideg-1),in)=tranmatrix(kp(in+1),kp(in+2),tai);
        tdown(5+8*(ideg-1):8+8*(ideg-1),5+8*(ideg-1):8+8*(ideg-1),in)=tranmatrix(kp(in+1),kp(in+2),tai);
    end
end

rup0 = zeros(4*ndeg,4*ndeg,nl);  %r10,r21,...rn,n-1,
rdown0 = zeros(4*ndeg,4*ndeg,nl);  %r12,r23,...rn,n+1
tup0 = zeros(4*ndeg,4*ndeg,nl);  %t10,t21,...tn,n-1,
tdown0 = zeros(4*ndeg,4*ndeg,nl);  %t12,t23,...tn,n+1
% kp=k0,k1,k2,...kn
for in=1:nl
    for ideg=1:ndeg
        tai=acos(mu(ideg));
        rup0(1+4*(ideg-1):4+4*(ideg-1),1+4*(ideg-1):4+4*(ideg-1),in)=refmatrix(kp(in+1),kp(in),tai);
        rdown0(1+4*(ideg-1):4+4*(ideg-1),1+4*(ideg-1):4+4*(ideg-1),in)=refmatrix(kp(in+1),kp(in+2),tai);
        
        tup0(1+4*(ideg-1):4+4*(ideg-1),1+4*(ideg-1):4+4*(ideg-1),in)=tranmatrix(kp(in+1),kp(in),tai);
        tdown0(1+4*(ideg-1):4+4*(ideg-1),1+4*(ideg-1):4+4*(ideg-1),in)=tranmatrix(kp(in+1),kp(in+2),tai);
    end
end
