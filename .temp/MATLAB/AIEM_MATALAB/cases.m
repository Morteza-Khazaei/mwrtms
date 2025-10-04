clear
clc
%%  backscattering case
theta_i=0:1:60;
phi_s=179.9999;
kl=2;
ks=0.4;
err=12;
eri=1.8;
itype=[1 2];
er=err-eri*1j;

for n=1:1:length(itype)
for i=1:1:length(theta_i)
theta_s=theta_i(i);
[sigHH, sigVV, sigHV ,sigVH]=AIEM(theta_i(i), theta_s,phi_s, kl, ks, err, eri, itype(n));
HH(n,i)=sigHH;
VV(n,i)=sigVV;
HV(n,i)=sigHV;
VH(n,i)=sigVH;
end
end

angle=0:1:60;
figure(1)
plot(angle,HH(1,:),'color','r','LineWidth',2)
hold on
plot(angle,VV(1,:),'color','r','LineWidth',2,'LineStyle','--')
plot(angle,HH(2,:),'color','b','LineWidth',2)
plot(angle,VV(2,:),'color','b','LineStyle','--','LineWidth',2)

legend({'HH_G_a_u_s_s','VV_G_a_u_s_s','HH_e_x_p'...
   'VV_e_x_p'},'Location','southwest','FontSize',12)
grid on

set(gca,'xtickLabel',{'0\circ','10\circ','20\circ','30\circ','40\circ','50\circ'...
    ,'60\circ'})
xlabel('Incident angle','Fontsize',16,'FontWeight','bold')
ylabel('Backscattering Coefficients (dB)','Fontsize',16,'FontWeight','bold')

%%  bistatic case
theta_i=30;
theta_s=60;
phi_s=0:1:180;
kl=2;
ks=0.4;
err=12;
eri=1.8;
itype=[1 2];
er=err-eri*1j;

for n=1:1:length(itype)
for i=1:1:length(phi_s)
[sigHH, sigVV, sigHV ,sigVH]=AIEM(theta_i, theta_s,phi_s(i), kl, ks, err, eri, itype(n));
HH(n,i)=sigHH;
VV(n,i)=sigVV;
HV(n,i)=sigHV;
VH(n,i)=sigVH;
end
end

angle=0:1:180;
figure(2)
plot(angle,HH(1,:),'color','r','LineWidth',2)
hold on
plot(angle,VV(1,:),'color','r','LineWidth',2,'LineStyle','--')
plot(angle,HH(2,:),'color','b','LineWidth',2)
plot(angle,VV(2,:),'color','b','LineStyle','--','LineWidth',2)

legend({'HH_G_a_u_s_s','VV_G_a_u_s_s','HH_e_x_p'...
   'VV_e_x_p'},'Location','southwest','FontSize',12)
grid on

set(gca,'xtickLabel',{'0\circ','20\circ','40\circ','60\circ','80\circ','100\circ'...
    ,'120\circ','140\circ','160\circ','180\circ'})
xlabel('Scattering Azimuth Angle','Fontsize',16,'FontWeight','bold')
ylabel('Bistatic Scattering Coefficients (dB)','Fontsize',16,'FontWeight','bold')
