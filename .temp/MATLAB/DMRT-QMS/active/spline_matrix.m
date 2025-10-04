function v= spline_matrix(x,xx)
%SPLINE Cubic spline data interpolation.
%   YY = SPLINE(X,Y,XX) uses cubic spline interpolation to find YY, the values
%   of the underlying function Y at the points in the vector XX.  The vector X
%   specifies the points at which the data Y is given.  If Y is a matrix, then
%   the data is taken to be vector-valued and interpolation is performed for
%   each column of Y and YY will be length(XX)-by-size(Y,2).
%
%   PP = SPLINE(X,Y) returns the piecewise polynomial form of the cubic spline
%   interpolant for later use with PPVAL and the spline utility UNMKPP.
%
%   Ordinarily, the not-a-knot end conditions are used. However, if Y contains
%   two more values than X has entries, then the first and last value in Y are
%   used as the endslopes for the cubic spline.  Namely:
%      f(X) = Y(:,2:end-1),   df(min(X)) = Y(:,1),   df(max(X)) = Y(:,end)
%
%   Example:
%   This generates a sine curve, then samples the spline over a finer mesh:
%       x = 0:10;  y = sin(x);
%       xx = 0:.25:10;
%       yy = spline(x,y,xx);
%       plot(x,y,'o',xx,yy)
%
%   Example:
%   This illustrates the use of clamped or complete spline interpolation where
%   end slopes are prescribed. Zero slopes at the ends of an interpolant to the
%   values of a certain distribution are enforced:
%      x = -4:4; y = [0 .15 1.12 2.36 2.36 1.46 .49 .06 0];
%      cs = spline(x,[0 y 0]);
%      xx = linspace(-4,4,101);
%      plot(x,y,'o',xx,ppval(cs,xx),'-');
%
%   See also INTERP1, PPVAL, SPLINES (The Spline Toolbox).

%   Carl de Boor 7-2-86
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 5.18 $  $Date: 2002/06/06 13:39:51 $

% Generate the cubic spline interpolant in ppform, depending on the
% number of data points (and usually using the not-a-knot end condition).

output=[];
n=length(x);
if n<2, error('There should be at least two data points.'), end

if any(diff(x)<0), [x,ind]=sort(x); else, ind=1:n; end

x=x(:); dx = diff(x);
if all(dx)==0, error('The data abscissae should be distinct.'), end

%[yd,yn] = size(y); % if Y happens to be a column matrix, change it to 
% the expected row matrix.
%if yn==1, yn=yd; y=reshape(y,1,yn); yd=1; end
% 
% if yn==n
notaknot = 1;
% elseif yn==n+2
%     notaknot = 0; endslopes = y(:,[1 n+2]).'; y(:,[1 n+2])=[];
% else
%     error('Abscissa and ordinate vector should be of the same length.')
% end

% yi=y(:,ind).'; dd = ones(1,yd);
dd = ones(1,1);
dx = diff(x); %divdif = diff(yi)./dx(:,dd);

dx1=diag(-1./dx); dx2=diag(1./dx); divy=[dx1,zeros(n-1,1)]+[zeros(n-1,1),dx2];
if n==2
    if notaknot, % the interpolant is a straight line
        pp=mkpp(x.',[divdif.' yi(1,:).'],yd);
    else         % the interpolant is the cubic Hermite polynomial
        divdif2 = diff([endslopes(1,:);divdif;endslopes(2,:)])./dx([1 1],dd);
        pp = mkpp(x,...
            [(diff(divdif2)./dx(1,dd)).' ([2 -1]*divdif2).' ...
                endslopes(1,:).' yi(1,:).'],yd);
    end
elseif n==3&notaknot, % the interpolant is a parabola
    yi(2:3,:)=divdif;
    yi(3,:)=diff(divdif)/(x(3)-x(1));
    yi(2,:)=yi(2,:)-yi(3,:)*dx(1);
    pp = mkpp([x(1),x(3)],yi([3 2 1],:).',yd);
else % set up the sparse, tridiagonal, linear system for the slopes at  X .
    %b=zeros(n,yd);
    %b(2:n-1,:)=3*(dx(2:n-1,dd).*divdif(1:n-2,:)+dx(1:n-2,dd).*divdif(2:n-1,:));
    if notaknot
        x31=x(3)-x(1);xn=x(n)-x(n-2);
        % b(1,:)=((dx(1)+2*x31)*dx(2)*divdif(1,:)+dx(1)^2*divdif(2,:))/x31;
        % b(n,:)=...
        % (dx(n-1)^2*divdif(n-2,:)+(2*xn+dx(n-1))*dx(n-2)*divdif(n-1,:))/xn;
    else
        x31 = 0; xn = 0; b([1 n],:) = dx([2 n-2],dd).*endslopes;
    end
    
    t=zeros(n,n);
    for id=2:n-1
        t(id,id-1:id+1)=[-3*dx(id)/dx(id-1),3*(dx(id)/dx(id-1)-dx(id-1)/dx(id)),3*dx(id-1)/dx(id)];
    end
    
    x31=x(3)-x(1);xn=x(n)-x(n-2);
    t(1,1:3)=[-(dx(1)+2*x31)*dx(2)/dx(1)/x31,((dx(1)+2*x31)*dx(2)/dx(1)-dx(1)^2/dx(2))/x31,dx(1)^2/x31/dx(2)];
    t(n,n-2:n)=[-dx(n-1)^2/dx(n-2)/xn,(dx(n-1)^2/dx(n-2)-(2*xn+dx(n-1))*dx(n-2)/dx(n-1))/xn,(2*xn+dx(n-1))*dx(n-2)/dx(n-1)/xn];
    
    
    c = spdiags([ [dx(2:n-1);xn;0] ...
            [dx(2);2*[dx(2:n-1)+dx(1:n-2)];dx(n-2)] ...
            [0;x31;dx(1:n-2)] ],[-1 0 1],n,n);
    
    % sparse linear equation solution for the slopes
    mmdflag = spparms('autommd');
    spparms('autommd',0);
    sy=c\t;
    %s=sy*y;
    spparms('autommd',mmdflag);
    % convert to pp form
    %c4=(s(1:n-1,:)+s(2:n,:)-2*divdif(1:n-1,:))./dx(:,dd);
    %c3=(divdif(1:n-1,:)-s(1:n-1,:))./dx(:,dd) - c4;
    cc4=(sy(1:n-1,:)+sy(2:n,:)-2*divy(1:n-1,:));
    cc3=(divy(1:n-1,:)-sy(1:n-1,:)) - cc4;
    
    c4y=zeros(n-1,n);
    c3y=zeros(n-1,n);
    for ic=1:n-1
        c4y(ic,:)=cc4(ic,:)/dx(ic);
        c3y(ic,:)=cc3(ic,:)/dx(ic);
        c41y(ic,:)=cc4(ic,:)/dx(ic)/dx(ic);
    end
    
    [mx,nx] = size(xx); lx = mx*nx; xs = reshape(xx,1,lx);  
    % if necessary, sort xx 
    tosort=0;
    if any(diff(xs)<0)
        tosort=1;[xs,ix]=sort(xs);
    end
    
    l=n-1;
    x = reshape(x,1,n); 
    [ignored,index] = sort([x(1:l) xs]);
    index = max([find(index>l)-(1:lx);ones(1,lx)]);
    
    % now go to local coordinates ...
    xs = xs-x(index);
    
    C1=c41y;
    C2=c3y;
    C3=sy(1:n-1,:);
    C4=[eye(n-1) zeros(n-1,1)];
    
    C1=C1(index,:);
    C2=C2(index,:);
    C3=C3(index,:);
    C4=C4(index,:);
    v=zeros(lx,n);
    
    for ix=1:lx
        v(ix,:)= C1(ix,:)*xs(ix)^3+C2(ix,:)*xs(ix)^2+C3(ix,:)*xs(ix)+C4(ix,:);
    end
    
end

np=length(x);n=length(xx);
[mv,nv]=size(v); lv=mv*nv; vs= reshape(v',1,lv);
indv=lv:-1:1;
vinv=vs(indv);
v=reshape(vinv,np,n);
v=v';
