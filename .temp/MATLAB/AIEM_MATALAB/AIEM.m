function [HH, VV, HV ,VH]=AIEM(theta_i, theta_s,phi_s, kl, ks, err, eri, itype)
%                
%   Computes the scattering from 3 types of correlated surface(Gauss&exponential&1.5-power)
%   by using Advanced Integral Equation Model (AIEM)

%   INPUT:
%   theta_i= incident angle in degree.
%   phi_s=scattering azimuth angle in deg while  incident azimuth angle
%               is 0 deg.
%   kl  = normalized surface correlation length multiplies by wave number k.
%   ks = normalized surface rms height multiplies by wave number k.
%   err = the real part of surface relative dielectric constant.
%   eri = the imaginary part of surface relative dielectric constant.
%   type= select what type of surface correlation
%   itype =1 Gaussian correlation function
%   itype =2 exponential correlation function
%   itype =3 transformed exponential correlation (1.5-power)
%
%   OUTPUT:
%
%   [HH, VV, HV ,VH]=scattering coefficient of HH-pol, VV-pol, HV-pol,
%   and VH-pol

    global si  sis cs css sfs csfs si2 sis2 cs2 css2 rv rh rhv ks2 er 
    phi_i=0;
    er=err+eri*1j;
    ur=1;
    ityp=itype;
    surface_1=ityp;
    theta_i=theta_i/180*pi;
    theta_s=theta_s/180*pi;
    phi_i=phi_i/180*pi;
    phi_s=phi_s/180*pi;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    si=sin(theta_i);
    sis=sin(theta_s);
    sfs = sin(phi_s);
    cs=cos(theta_i);
    css=cos(theta_s);  
    csfs = cos(phi_s);
    
    cs2 = cs.*cs;
    css2 = css.*css;
    si2 = si.*si;
    sis2 = sis.*sis;
    ks2= ks.*ks;
    kl2 = kl.*kl;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% compute roughness spectrum spectra_1(n) 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    torlant = 1.0e-16;
    iterm = 1;
    tempold = 0.0;
    temp =(ks2.*(cs+css).^2.0);
    while(abs(temp-tempold)>torlant) 
         tempold = temp;
         iterm = iterm + 1;
         fiterm =iterm;
         temp = tempold.*(ks2.*(cs+css).^2.0)./fiterm;
    end
    
    for n = 1: 1000
         spectra_1(n) = 0.0;
    end

    for n = 1: iterm
         fn = n;
         K = kl.*sqrt((sis.*csfs-si*cos(phi_i)).^2.0+(sis.*sfs-si*sin(phi_i)).^2.0);%指定空间频率成分
         itype=num2str(surface_1);
         switch lower(itype)
              case {'1'}
%------ Gaussian correlated surface -----------c
                   spectra_1(1,n) = kl2.*exp(-K.*K./(4.0.*fn))./(2.*fn);
             case {'2'}
%------ exponential correlated surface -----------c
                   spectra_1(1,n) =(kl./fn).^2.0.*(1.0+(K./fn).^2.0).^(-1.5);
            case {'3'}
 %------ 1.5 power surface -----------c      
                  e = 1.5.*fn - 1.0 ;
                  y = 1.5.*fn;   
                  gam =log(gamma(y));%gamma function (1.5n)
                  if K==0.0
                     spectra_1(n) = kl.*kl./(3.0.*fn-2.0);
                  else
                          m = 1.5.*fn -1.0;
                          bk = log(besselk(-m,K));
                          out = kl.*kl.*(K./2.0 ).^e; 
                          spectra_1(n) = out.*exp(bk-gam);
                  end
         end
    end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
% reflection coefficients based on the incident angle  
% ==> R(theta_i)                                          
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    stem = sqrt(er.*ur-si2);
    rvi=(er.*cs-stem)./(er.*cs+stem);
    rhi=(ur.*cs-stem)./(ur.*cs+stem);
    rvhi=(rvi-rhi)./2.0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% reflection coefficients based on the specular angle  
% ==> R(theta_specular)                                 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    csl = sqrt(1.0 +cs.*css-si.*sis.*csfs)./sqrt(2.0);
    sil = sqrt(1.0 -csl.*csl);
    steml = sqrt(er.*ur-sil.*sil);
    rvl =(er.*csl-steml)./(er.*csl+steml);
    rhl =(ur.*csl-steml)./(ur.*csl+steml);
    rvhl =(rvl-rhl)./2.0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% reflection coefficients based on the transition function  
% ==> R_transition (T.D.Wu&A.K.fung)                       
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rv0 =(sqrt(er)-1.0)./(sqrt(er)+1.0);
    rh0 =-(sqrt(er)-1.0)./(sqrt(er)+1.0);
    Ftv = 8.0.*(rv0.^2).*si2.*(cs+sqrt(er-si2))/(cs.*(sqrt(er-si2)));
    Fth = -8.0.*(rh0.^2).*si2.*(cs+sqrt(er-si2))/(cs.*(sqrt(er-si2)));
    st0v =1.0./((abs(1.0 +8.0 .*rv0./(cs.*Ftv))).^2.0);
    st0h =1.0./((abs(1.0 +8.0 .*rv0./(cs.*Fth))).^2.0);
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    temp1 = 1.0;

    for n = 1: iterm
        fn =n;
        temp1= temp1.*(1./fn);
        sum1= sum1 + temp1.*((ks.*cs).^(2.0 .*fn)).*spectra_1(n);
        sum2= sum2 + temp1.*((ks.*cs).^(2.0 .*fn))*(abs(Ftv+2.0 .^(fn+2.0 ).*rv0./cs./(exp((ks.*cs).^2.0 ))).^2.0 ).*spectra_1(n);
        sum3= sum3 + temp1.*((ks.*cs).^(2.0 .*fn))*(abs(Fth+2.0 .^(fn+2.0 ).*rv0./cs.*(exp(-(ks.*cs).^2.0 ))).^2.0 ).*spectra_1(n);
    end

    stv=(abs(Ftv).^2.0 ).*sum1./sum2;
    sth=(abs(Fth).^2.0 ).*sum1./sum3;
    tfv = 1.0-stv./st0v;
    tfh = 1.0-sth./st0h;
    if (tfv<0.0 )
        tfv = 0.0;
    end
    if(tfh<0.0 )
        tfh = 0.0;
    end
    rvtran = rvi +(rvl-rvi).*tfv;
    rhtran = rhi +(rhl-rhi).*tfh;
    rvhtran =(rvtran-rhtran)./2.0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     Reflection coefficients rv, rh, rhv  for kirchhoff field coefficients            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rh =rhtran;
    rv = rvtran;
    rhv = rvhtran;
%----------------------------------------------------
%     kirchhoff field coefficients fvv fhh fhv fvh 
%----------------------------------------------------
    zxx = -(sis.*csfs-si)./(css+cs);
    zyy = -(sis.*sfs)./(css+cs);
    d2 = sqrt((zxx.*cs-si).^2.0 +zyy.^2.0 );
    hsnv = -(cs.*csfs+si.*(zxx.*csfs+zyy.*sfs));
    vsnh = css.*csfs - zxx.*sis;
    hsnh = -sfs;
    vsnv = zyy.*cs.*sis + css.*(zyy.*csfs.*si-(cs+zxx.*si).*sfs);
    hsnt =(-(cs2+si2).*sfs.*(-si+cs.*zxx)+csfs.*(cs+si.*zxx).*zyy+si.*sfs.*(zyy.^2))./d2;
    hsnd =(-(cs+si.*zxx).*(-csfs.*si+cs.*csfs.*zxx+cs.*sfs.*zyy))./d2;
    vsnt =((cs2+si2).*(-si+cs.*zxx).*(csfs.*css-sis.*zxx)+css.*sfs.*(cs+si.*zxx).*zyy-(csfs.*css.*si+cs.*sis).*(zyy.^2))./d2;
    vsnd = -(cs+si.*zxx).*(si.*sis.*zyy-css.*(si.*sfs-cs.*sfs.*zxx+cs.*csfs.*zyy))./d2;

    fhh =(1.0 -rh).*hsnv +(1.0 +rh).*vsnh -(hsnt+vsnd).*(rh+rv).*(zyy./d2);
    fvv = -((1.0 -rv).*hsnv+(1.0 +rv).*vsnh) +(hsnt+vsnd).*(rh+rv).*(zyy./d2);
    fhv = -(1.0 +rv).*hsnh +(1.0 -rv).*vsnv +(hsnd-vsnt).*(rh+rv).*(zyy./d2);
    fvh = -(1.0 +rh).*hsnh +(1.0 -rh).*vsnv +(hsnd-vsnt).*(rh+rv).*(zyy./d2);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    compute kirchhoff term                         
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sum = 0.0;
    temp = 1.0;
    for n = 1: iterm

        fn = n;
        temp = temp.*(ks2.*(cs+css).^2.0)./fn;
        sum = sum + temp.*spectra_1(n);
    end
    expk = exp(-ks2.*(css+cs).^2.0).*sum;
    kterm(1) = 0.5.*expk.*abs(fvv).^2;
    kterm(2) = 0.5.*expk.*abs(fhh).^2;
    kterm(3) = 0.5.*expk.*abs(fhv).^2;
    kterm(4) = 0.5.*expk.*abs(fvh).^2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     Reflection coefficients rv, rh, rhv for complementary field coefficients       
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rh = rhi;
    rv = rvi;
    rhv = rvhi;
    qq = cs;
    qqt = sqrt(er-si2);
    qqs = css;
    qqts = sqrt(er-sis2);

    qq1 = qq;
    qq2 = qqs;
    qq3 = qqt;
    qq4 = qqts;
    qq5 = qqt;
    qq6 = qqts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%!!!!!Fvv   !!!!!!!!!
    Fvaupi = favv(-si, 0.0 , qq1, qq1, qq).*expal(qq1);%qq1向上
    Fvadni = favv(-si, 0.0 , -qq1, -qq1, qq).*expal(-qq1);%-qq1向下
    Fvaups = favv(-sis.*csfs, -sis.*sfs, qq2, qq2, qqs).*expal(qq2);
    Fvadns = favv(-sis.*csfs, -sis.*sfs, -qq2, -qq2, qqs).*expal(-qq2);
    Fvbupi = fbvv(-si, 0.0 , qq3, qq5, qqt).*expal(qq5);
    Fvbdni = fbvv(-si, 0.0 , -qq3, -qq5, qqt).*expal(-qq5);
    Fvbups = fbvv(-sis.*csfs, -sis.*sfs, qq4, qq6, qqts).*expal(qq6);
    Fvbdns = fbvv(-sis.*csfs, -sis.*sfs, -qq4, -qq6, qqts).*expal(-qq6);
%!!!!Fhh!!!!!!!!!
    Fhaupi = fahh(-si, 0.0 , qq1, qq1, qq).*expal(qq1);
    Fhadni = fahh(-si, 0.0 , -qq1, -qq1, qq).*expal(-qq1);
    Fhaups = fahh(-sis.*csfs, -sis.*sfs, qq2, qq2, qqs).*expal(qq2);
    Fhadns = fahh(-sis.*csfs, -sis.*sfs, -qq2, -qq2, qqs).*expal(-qq2);
    Fhbupi = fbhh(-si, 0.0 , qq3, qq5, qqt).*expal(qq5);
    Fhbdni = fbhh(-si, 0.0 , -qq3, -qq5, qqt).*expal(-qq5);
    Fhbups = fbhh(-sis.*csfs, -sis.*sfs, qq4, qq6, qqts).*expal(qq6);
    Fhbdns = fbhh(-sis.*csfs, -sis.*sfs, -qq4, -qq6, qqts).*expal(-qq6);
%!!!!Fhv!!!!!!!!!
    Fhvaupi = fahv(-si, 0.0 , qq1, qq1, qq).*expal(qq1);
    Fhvadni = fahv(-si, 0.0 , -qq1, -qq1, qq).*expal(-qq1);
    Fhvaups = fahv(-sis.*csfs, -sis.*sfs, qq2, qq2, qqs).*expal(qq2);
    Fhvadns = fahv(-sis.*csfs, -sis.*sfs, -qq2, -qq2, qqs).*expal(-qq2);
    Fhvbupi = fbhv(-si, 0.0 , qq3, qq5, qqt).*expal(qq5);
    Fhvbdni = fbhv(-si, 0.0 , -qq3, -qq5, qqt).*expal(-qq5);
    Fhvbups = fbhv(-sis.*csfs, -sis.*sfs, qq4, qq6, qqts).*expal(qq6);
    Fhvbdns = fbhv(-sis.*csfs, -sis.*sfs, -qq4, -qq6, qqts).*expal(-qq6);
%!!!!Fvh!!!!!!!!!
    Fvhaupi = favh(-si, 0.0 , qq1, qq1, qq).*expal(qq1);
    Fvhadni = favh(-si, 0.0 , -qq1, -qq1, qq).*expal(-qq1);
    Fvhaups = favh(-sis.*csfs, -sis.*sfs, qq2, qq2, qqs).*expal(qq2);
    Fvhadns = favh(-sis.*csfs, -sis.*sfs, -qq2, -qq2, qqs).*expal(-qq2);
    Fvhbupi = fbvh(-si, 0.0 , qq3, qq5, qqt).*expal(qq5);
    Fvhbdni = fbvh(-si, 0.0 , -qq3, -qq5, qqt).*expal(-qq5);
    Fvhbups = fbvh(-sis.*csfs, -sis.*sfs, qq4, qq6, qqts).*expal(qq6);
    Fvhbdns = fbvh(-sis.*csfs, -sis.*sfs, -qq4, -qq6, qqts).*expal(-qq6);

%--------------------------------------------------------------c
%    Compute scattering coefficients!
%--------------------------------------------------------------    
    for n = 1: iterm

        fn = n;
        Ivv(n) =((cs+css).^fn).*fvv.*exp(-ks2.*cs.*css)+(0.25 ).*... 
            (Fvaupi.*((css-qq1).^fn)+Fvadni.*((css+qq1).^fn)...
            +Fvaups.*((cs+qq2).^fn)+Fvadns.*((cs-qq2).^fn)...
            +Fvbupi.*((css-qq5).^fn)+Fvbdni.*((css+qq5).^fn)...
            +Fvbups.*((cs+qq6).^fn)+Fvbdns.*((cs-qq6).^fn));

        Ihh(n) =((cs+css).^fn).*fhh.*exp(-ks2.*cs.*css) +(0.25 ).*...
            (Fhaupi.*((css-qq1).^fn)+Fhadni.*((css+qq1).^fn)...
            +Fhaups.*((cs+qq2).^fn)+Fhadns.*((cs-qq2).^fn)...
            +Fhbupi.*((css-qq5).^fn)+Fhbdni.*((css+qq5).^fn)...
            +Fhbups.*((cs+qq6).^fn)+Fhbdns.*((cs-qq6).^fn));

        Ihv(n) =((cs+css).^fn).*fhv.*exp(-ks2.*cs.*css) +(0.25 ).*...
            (Fhvaupi.*((css-qq1).^fn)+Fhvadni.*((css+qq1).^fn)...
            +Fhvaups.*((cs+qq2).^fn)+Fhvadns.*((cs-qq2).^fn)...
            +Fhvbupi.*((css-qq5).^fn)+Fhvbdni.*((css+qq5).^fn)...
            +Fhvbups.*((cs+qq6).^fn)+Fhvbdns.*((cs-qq6).^fn));

        Ivh(n) =((cs+css).^fn).*fvh.*exp(-ks2.*cs.*css) +(0.25 ).*...
            (Fvhaupi.*((css-qq1).^fn)+Fvhadni.*((css+qq1).^fn)...
            +Fvhaups.*((cs+qq2).^fn)+Fvhadns.*((cs-qq2).^fn)...
            +Fvhbupi.*((css-qq5).^fn)+Fvhbdni.*((css+qq5).^fn)...
            +Fvhbups.*((cs+qq6).^fn)+Fvhbdns.*((cs-qq6).^fn));

        CIvv(n) = conj(Ivv(n));
        CIhh(n) = conj(Ihh(n));
        CIhv(n) = conj(Ihv(n));
        CIvh(n) = conj(Ivh(n));

    end

    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    temp = 1.0;
    for n = 1: iterm

        fn = n;
        temp = temp.*(ks2./fn);
        sum1 = sum1 + temp.*(Ivv(n).*CIvv(n)).*spectra_1(n);
        sum2 = sum2 + temp.*(Ihh(n).*CIhh(n)).*spectra_1(n);
        sum3 = sum3 + temp.*(Ihv(n).*CIhv(n)).*spectra_1(n);
        sum4 = sum4 + temp.*(Ivh(n).*CIvh(n)).*spectra_1(n);
    end
    allterm(1) =(0.5 ).*exp(-ks2.*(cs2+css2)).*sum1;
    allterm(2) =(0.5 ).*exp(-ks2.*(cs2+css2)).*sum2;
    allterm(3) =(0.5 ).*exp(-ks2.*(cs2+css2)).*sum3;
    allterm(4) =(0.5 ).*exp(-ks2.*(cs2+css2)).*sum4;
    sigma0(1) =real(allterm(1));
    sigma0(2) = real(allterm(2));
    sigma0(3) = real(allterm(3));
    sigma0(4) = real(allterm(4));
    VV=10*log10(sigma0(1));
    HH=10*log10(sigma0(2));
    HV=10*log10(sigma0(3));
    VH=10*log10(sigma0(4));
end
%end of main function program


function [expalresult] =expal(q)
    global  cs css ks2;
    expalresult = exp(-ks2.*(q.^2.0 -q.*(css-cs)));%
end 


function [fahhresult]=fahh(u, v, q, qslp, qfix)
    global si  sis cs css sfs csfs rh;
    kxu = si + u;
    ksxu = sis.*csfs + u;
    kyv = v;
    ksyv = sis.*sfs + v;
    if(abs(real(css-qslp))<0.0000000001)
        zx = 0.0 ;
        zy = 0.0 ;
    else
        zx =(-ksxu)./(css-qslp);
        zy = -(ksyv)./(css-qslp);
    end
    if(abs(real(cs+qslp))<0.0000000001)
        zxp = 0.0 ;
        zyp = 0.0 ;
    else
        zxp =(kxu)./(cs+qslp);
        zyp =(kyv)./(cs+qslp);
    end
    c1 = -csfs.*(-1.0-zx.*zxp) + sfs.*zxp.*zy;
    c2 = -csfs.*(-cs.*q-cs.*u.*zx-q.*si.*zxp-si.*u.*zx.*zxp-cs.*v.*zyp-si.*v.*zx.*zyp) ...
        + sfs.*(cs.*u.*zy+si.*u.*zxp.*zy+q.*si.*zyp-cs.*u.*zyp+si.*v.*zy.*zyp);
    c3 = -csfs.*(si.*u-q.*si.*zx-cs.*u.*zxp+cs.*q.*zx.*zxp) ...
        + sfs.*(-si.*v+cs.*v.*zxp+q.*si.*zy-cs.*q.*zxp.*zy);
    c4 = -css.*sfs.*(-si.*zyp+cs.*zx.*zyp) - csfs.*css.*(-cs-si.*zxp-cs.*zy.*zyp) + sis.*(-cs.*zx-si.*zx.*zxp-si.*zy.*zyp);
    c5 = -css.*sfs.*(-v.*zx+v.*zxp) - csfs.*css.*(q+u.*zxp+v.*zy) + sis.*(q.*zx+u.*zx.*zxp+v.*zxp.*zy);
    c6 = -css.*sfs.*(-u.*zyp+q.*zx.*zyp) - csfs.*css.*(v.*zyp-q.*zy.*zyp) + sis.*(v.*zx.*zyp-u.*zy.*zyp);

%
    rph = 1.0 + rh;
    rmh = 1.0 - rh;
    ah = rph./qfix;
    bh = rmh./qfix;
    fahhresult = -bh.*(-rph.*c1+rmh.*c2+rph.*c3) - ah.*(rmh.*c4+rph.*c5+rmh.*c6);
end %function fahh



function [fahvresult]=fahv(u, v, q, qslp, qfix)
      global si  sis cs css sfs csfs rhv;
    kxu = si + u;
    ksxu = sis.*csfs + u;
    kyv = v;
    ksyv = sis.*sfs + v;
    if(abs(real(css-qslp))<0.0000000001)
        zx = 0.0 ;
        zy = 0.0 ;
    else
        zx =(-ksxu)./(css-qslp);
        zy = -(ksyv)./(css-qslp);
    end
    if(abs(real(cs+qslp))<0.0000000001)
        zxp = 0.0 ;
        zyp = 0.0 ;
    else
        zxp =(kxu)./(cs+qslp);
        zyp =(kyv)./(cs+qslp);
    end

    b1 = -css.*sfs.*(-1.0-zx.*zxp) - sis.*zy - csfs.*css.*zxp.*zy;
    b2 = -css.*sfs.*(-cs.*q-cs.*u.*zx-q.*si.*zxp-si.*u.*zx.*zxp-cs.*v.*zyp-si.*v.*zx.*zyp) ...
        + sis.*(-cs.*q.*zy-q.*si.*zxp.*zy+q.*si.*zx.*zyp-cs.*u.*zx.*zyp-cs.*v.*zy.*zyp) ...
        - csfs.*css.*(cs.*u.*zy+si.*u.*zxp.*zy+q.*si.*zyp-cs.*u.*zyp+si.*v.*zy.*zyp);
    b3 = -css.*sfs.*(si.*u-q.*si.*zx-cs.*u.*zxp+cs.*q.*zx.*zxp) - ...
        csfs.*css.*(-si.*v+cs.*v.*zxp+q.*si.*zy-cs.*q.*zxp.*zy) ...
        + sis.*(-si.*v.*zx+cs.*v.*zx.*zxp+si.*u.*zy-cs.*u.*zxp.*zy);
    b4 = -csfs.*(-si.*zyp+cs.*zx.*zyp) + sfs.*(-cs-si.*zxp-cs.*zy.*zyp);
    b5 = -csfs.*(-v.*zx+v.*zxp) + sfs.*(q+u.*zxp+v.*zy);
    b6 = -csfs.*(-u.*zyp+q.*zx.*zyp) + sfs.*(v.*zyp-q.*zy.*zyp);

    rp = 1.0 + rhv;
    rm = 1.0 - rhv;
    a = rp./qfix;
    b = rm./qfix;
    fahvresult = b.*(rp.*b1-rm.*b2-rp.*b3) + a.*(rm.*b4+rp.*b5+rm.*b6);
end %function fahv


function [favhresult]=favh(u, v, q, qslp, qfix)
    global si  sis cs css sfs csfs rhv;
    kxu = si + u;
    ksxu = sis.*csfs + u;
    kyv = v;
    ksyv = sis.*sfs + v;
    if(abs(real(css-qslp))<0.0000000001)
        zx = 0.0 ;
        zy = 0.0 ;
    else
        zx =(-ksxu)./(css-qslp);
        zy = -(ksyv)./(css-qslp);
    end
    if(abs(real(cs+qslp))<0.0000000001)
        zxp = 0.0 ;
        zyp = 0.0 ;
    else
        zxp =(kxu)./(cs+qslp);
        zyp =(kyv)./(cs+qslp);
    end

    b1 = -css.*sfs.*(-1.0-zx.*zxp) - sis.*zy - csfs.*css.*zxp.*zy;
    b2 = -css.*sfs.*(-cs.*q-cs.*u.*zx-q.*si.*zxp-si.*u.*zx.*zxp-cs.*v.*zyp-si.*v.*zx.*zyp) ...
        + sis.*(-cs.*q.*zy-q.*si.*zxp.*zy+q.*si.*zx.*zyp-cs.*u.*zx.*zyp-cs.*v.*zy.*zyp)...
        - csfs.*css.*(cs.*u.*zy+si.*u.*zxp.*zy+q.*si.*zyp-cs.*u.*zyp+si.*v.*zy.*zyp);
    b3 = -css.*sfs.*(si.*u-q.*si.*zx-cs.*u.*zxp+cs.*q.*zx.*zxp)...
        - csfs.*css.*(-si.*v+cs.*v.*zxp+q.*si.*zy-cs.*q.*zxp.*zy) ...
        + sis.*(-si.*v.*zx+cs.*v.*zx.*zxp+si.*u.*zy-cs.*u.*zxp.*zy);
    b4 = -csfs.*(-si.*zyp+cs.*zx.*zyp) + sfs.*(-cs-si.*zxp-cs.*zy.*zyp);
    b5 = -csfs.*(-v.*zx+v.*zxp) + sfs.*(q+u.*zxp+v.*zy);
    b6 = -csfs.*(-u.*zyp+q.*zx.*zyp) + sfs.*(v.*zyp-q.*zy.*zyp);

    rp = 1.0 + rhv;
    rm = 1.0 - rhv;
    a = rp./qfix;
    b = rm./qfix;
    favhresult = b.*(rp.*b4+rm.*b5+rp.*b6) - a.*(-rm.*b1+rp.*b2+rm.*b3);
end %function favh

function [favvresult]=favv(u, v, q, qslp, qfix)
    global si  sis cs css sfs csfs rv;
    kxu = si + u;
    ksxu = sis.*csfs + u;
    kyv = v;
    ksyv = sis.*sfs + v;
    if(abs(real(css-qslp))<0.0000000001)
        zx = 0.0 ;
        zy = 0.0 ;
    else
        zx =(-ksxu)./(css-qslp);
        zy = -(ksyv)./(css-qslp);
    end
    if(abs(real(cs+qslp))<0.0000000001)
        zxp = 0.0 ;
        zyp = 0.0 ;
    else
        zxp =(kxu)./(cs+qslp);
        zyp =(kyv)./(cs+qslp);
    end
    c1 = -csfs.*(-1.0-zx.*zxp) + sfs.*zxp.*zy;
    c2 = -csfs.*(-cs.*q-cs.*u.*zx-q.*si.*zxp-si.*u.*zx.*zxp-cs.*v.*zyp-si.*v.*zx.*zyp) ...
        + sfs.*(cs.*u.*zy+si.*u.*zxp.*zy+q.*si.*zyp-cs.*u.*zyp+si.*v.*zy.*zyp);
    c3 = -csfs.*(si.*u-q.*si.*zx-cs.*u.*zxp+cs.*q.*zx.*zxp) ...
        + sfs.*(-si.*v+cs.*v.*zxp+q.*si.*zy-cs.*q.*zxp.*zy);
    c4 = -css.*sfs.*(-si.*zyp+cs.*zx.*zyp) - csfs.*css.*(-cs-si.*zxp-cs.*zy.*zyp) ...
        + sis.*(-cs.*zx-si.*zx.*zxp-si.*zy.*zyp);
    c5 = -css.*sfs.*(-v.*zx+v.*zxp) - csfs.*css.*(q+u.*zxp+v.*zy) ...
        + sis.*(q.*zx+u.*zx.*zxp+v.*zxp.*zy);
    c6 = -css.*sfs.*(-u.*zyp+q.*zx.*zyp) - csfs.*css.*(v.*zyp-q.*zy.*zyp)...
        + sis.*(v.*zx.*zyp-u.*zy.*zyp);
    rpv = 1.0 + rv;
    rmv = 1.0 - rv;
    av = rpv./qfix;
    bv = rmv./qfix;
    favvresult = bv.*(-rpv.*c1+rmv.*c2+rpv.*c3) + av.*(rmv.*c4+rpv.*c5+rmv.*c6);
end

function [fbhhresult]=fbhh(u, v, q, qslp, qfix)
    global si  sis cs css sfs csfs rh er;
    kxu = si + u;
    ksxu = sis.*csfs + u;
    kyv = v;
    ksyv = sis.*sfs + v;
    if(abs(real(css-qslp))<0.0000000001)
        zx = 0.0 ;
        zy = 0.0 ;
    else
        zx =(-ksxu)./(css-qslp);
        zy = -(ksyv)./(css-qslp);
    end
    if(abs(real(cs+qslp))<0.0000000001)
        zxp = 0.0 ;
        zyp = 0.0 ;
    else
        zxp =(kxu)./(cs+qslp);
        zyp =(kyv)./(cs+qslp);
    end
    c1 = -csfs.*(-1.0-zx.*zxp) + sfs.*zxp.*zy;
    c2 = -csfs.*(-cs.*q-cs.*u.*zx-q.*si.*zxp-si.*u.*zx.*zxp-cs.*v.*zyp-si.*v.*zx.*zyp) ...
        + sfs.*(cs.*u.*zy+si.*u.*zxp.*zy+q.*si.*zyp-cs.*u.*zyp+si.*v.*zy.*zyp);
    c3 = -csfs.*(si.*u-q.*si.*zx-cs.*u.*zxp+cs.*q.*zx.*zxp) ...
        + sfs.*(-si.*v+cs.*v.*zxp+q.*si.*zy-cs.*q.*zxp.*zy);
    c4 = -css.*sfs.*(-si.*zyp+cs.*zx.*zyp) - csfs.*css.*(-cs-si.*zxp-cs.*zy.*zyp)...
        + sis.*(-cs.*zx-si.*zx.*zxp-si.*zy.*zyp);
    c5 = -css.*sfs.*(-v.*zx+v.*zxp) - csfs.*css.*(q+u.*zxp+v.*zy) ...
        + sis.*(q.*zx+u.*zx.*zxp+v.*zxp.*zy);
    c6 = -css.*sfs.*(-u.*zyp+q.*zx.*zyp) - csfs.*css.*(v.*zyp-q.*zy.*zyp)...
        + sis.*(v.*zx.*zyp-u.*zy.*zyp);
    rph = 1.0 + rh;
    rmh = 1.0 - rh;
    ah = rph./qfix;
    bh = rmh./qfix;
    fbhhresult = ah.*(-rph.*c1.*er+rmh.*c2+rph.*c3) + bh.*(rmh.*c4+rph.*c5+rmh.*c6./er);
end %function fbhh


function [fbhvresult]=fbhv(u, v, q, qslp, qfix)
    global si  sis cs css sfs csfs rhv er
    kxu = si + u;
    ksxu = sis.*csfs + u;
    kyv = v;
    ksyv = sis.*sfs + v;
    if(abs(real(css-qslp))<0.0000000001)
        zx = 0.0 ;
        zy = 0.0 ;
    else
        zx =(-ksxu)./(css-qslp);
        zy = -(ksyv)./(css-qslp);
    end
    if(abs(real(cs+qslp))<0.0000000001)
        zxp = 0.0 ;
        zyp = 0.0 ;
    else
        zxp =(kxu)./(cs+qslp);
        zyp =(kyv)./(cs+qslp);
    end
    b1 = -css.*sfs.*(-1.0-zx.*zxp) - sis.*zy - csfs.*css.*zxp.*zy;
    b2 = -css.*sfs.*(-cs.*q-cs.*u.*zx-q.*si.*zxp-si.*u.*zx.*zxp-cs.*v.*zyp-si.*v.*zx.*zyp) ...
        + sis.*(-cs.*q.*zy-q.*si.*zxp.*zy+q.*si.*zx.*zyp-cs.*u.*zx.*zyp-cs.*v.*zy.*zyp) ...
        - csfs.*css.*(cs.*u.*zy+si.*u.*zxp.*zy+q.*si.*zyp-cs.*u.*zyp+si.*v.*zy.*zyp);
    b3 = -css.*sfs.*(si.*u-q.*si.*zx-cs.*u.*zxp+cs.*q.*zx.*zxp) ...
        - csfs.*css.*(-si.*v+cs.*v.*zxp+q.*si.*zy-cs.*q.*zxp.*zy) ...
        + sis.*(-si.*v.*zx+cs.*v.*zx.*zxp+si.*u.*zy-cs.*u.*zxp.*zy);
    b4 = -csfs.*(-si.*zyp+cs.*zx.*zyp) + sfs.*(-cs-si.*zxp-cs.*zy.*zyp);
    b5 = -csfs.*(-v.*zx+v.*zxp) + sfs.*(q+u.*zxp+v.*zy);
    b6 = -csfs.*(-u.*zyp+q.*zx.*zyp) + sfs.*(v.*zyp-q.*zy.*zyp);
    rp = 1.0 + rhv;
    rm = 1.0 - rhv;
    a = rp./qfix;
    b = rm./qfix;
    fbhvresult = a.*(-rp.*b1+rm.*b2+rp.*b3./er) - b.*(rm.*b4.*er+rp.*b5+rm.*b6);
end %function fbhv

function [fbvhresult]=fbvh(u, v, q, qslp, qfix)
    global si  sis cs css sfs csfs rhv er;
    kxu = si + u;
    ksxu = sis.*csfs + u;
    kyv = v;
    ksyv = sis.*sfs + v;
    if(abs(real(css-qslp))<0.0000000001)
        zx = 0.0 ;
        zy = 0.0 ;
    else
        zx =(-ksxu)./(css-qslp);
        zy = -(ksyv)./(css-qslp);
    end
    if(abs(real(cs+qslp))<0.0000000001)
        zxp = 0.0 ;
        zyp = 0.0 ;
    else
        zxp =(kxu)./(cs+qslp);
        zyp =(kyv)./(cs+qslp);
    end
    b1 = -css.*sfs.*(-1.0-zx.*zxp) - sis.*zy - csfs.*css.*zxp.*zy;
    b2 = -css.*sfs.*(-cs.*q-cs.*u.*zx-q.*si.*zxp-si.*u.*zx.*zxp-cs.*v.*zyp-si.*v.*zx.*zyp)...
        + sis.*(-cs.*q.*zy-q.*si.*zxp.*zy+q.*si.*zx.*zyp-cs.*u.*zx.*zyp-cs.*v.*zy.*zyp) ...
        - csfs.*css.*(cs.*u.*zy+si.*u.*zxp.*zy+q.*si.*zyp-cs.*u.*zyp+si.*v.*zy.*zyp);
    b3 = -css.*sfs.*(si.*u-q.*si.*zx-cs.*u.*zxp+cs.*q.*zx.*zxp) ...
        - csfs.*css.*(-si.*v+cs.*v.*zxp+q.*si.*zy-cs.*q.*zxp.*zy)...
        + sis.*(-si.*v.*zx+cs.*v.*zx.*zxp+si.*u.*zy-cs.*u.*zxp.*zy);
    b4 = -csfs.*(-si.*zyp+cs.*zx.*zyp) + sfs.*(-cs-si.*zxp-cs.*zy.*zyp);
    b5 = -csfs.*(-v.*zx+v.*zxp) + sfs.*(q+u.*zxp+v.*zy);
    b6 = -csfs.*(-u.*zyp+q.*zx.*zyp) + sfs.*(v.*zyp-q.*zy.*zyp);
    rp = 1.0 + rhv;
    rm = 1.0 - rhv;
    a = rp./qfix;
    b = rm./qfix;
    fbvhresult = -a.*(rp.*b4+rm.*b5+rp.*b6./er) + b.*(-rm.*b1.*er+rp.*b2+rm.*b3);
end %function fbvh

function [fbvvresult]=fbvv(u, v, q, qslp, qfix)
    global si  sis cs css sfs csfs rv er;
    kxu = si + u;
    ksxu = sis.*csfs + u;
    kyv = v;
    ksyv = sis.*sfs + v;
    if(abs(real(css-qslp))<0.0000000001)
        zx = 0.0 ;
        zy = 0.0 ;
    else
        zx =(-ksxu)./(css-qslp);
        zy = -(ksyv)./(css-qslp);
    end
    if(abs(real(cs+qslp))<0.0000000001)
        zxp = 0.0 ;
        zyp = 0.0 ;
    else
        zxp =(kxu)./(cs+qslp);
        zyp =(kyv)./(cs+qslp);
    end
    c1 = -csfs.*(-1.0-zx.*zxp) + sfs.*zxp.*zy;
    c2 = -csfs.*(-cs.*q-cs.*u.*zx-q.*si.*zxp-si.*u.*zx.*zxp-cs.*v.*zyp-si.*v.*zx.*zyp) ...
        + sfs.*(cs.*u.*zy+si.*u.*zxp.*zy+q.*si.*zyp-cs.*u.*zyp+si.*v.*zy.*zyp);
    c3 = -csfs.*(si.*u-q.*si.*zx-cs.*u.*zxp+cs.*q.*zx.*zxp) ...
        + sfs.*(-si.*v+cs.*v.*zxp+q.*si.*zy-cs.*q.*zxp.*zy);
    c4 = -css.*sfs.*(-si.*zyp+cs.*zx.*zyp) - csfs.*css.*(-cs-si.*zxp-cs.*zy.*zyp) ...
        + sis.*(-cs.*zx-si.*zx.*zxp-si.*zy.*zyp);
    c5 = -css.*sfs.*(-v.*zx+v.*zxp) - csfs.*css.*(q+u.*zxp+v.*zy) ...
        + sis.*(q.*zx+u.*zx.*zxp+v.*zxp.*zy);
    c6 = -css.*sfs.*(-u.*zyp+q.*zx.*zyp) - csfs.*css.*(v.*zyp-q.*zy.*zyp) ...
        + sis.*(v.*zx.*zyp-u.*zy.*zyp);

    rpv = 1.0 + rv;
    rmv = 1.0 - rv;
    av = rpv./qfix;
    bv = rmv./qfix;
    fbvvresult = av.*(rpv.*c1-rmv.*c2-rpv.*c3./er) - bv.*(rmv.*c4.*er+rpv.*c5+rmv.*c6);
end %function fbvv



