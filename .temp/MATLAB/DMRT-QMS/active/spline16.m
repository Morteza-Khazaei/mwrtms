
function v=spline16(epsilon_i,epsilon_t,mu)
% from epsilon_i to epsilon_t

n=length(mu);

if epsilon_i>epsilon_t
    % to a less dense media
    critical=cos(asin(sqrt(epsilon_t/epsilon_i)));  
    np=length(find(mu>critical));
    mup=zeros(1,np);
mup=sqrt(1-epsilon_i/epsilon_t*(1-mu(1:np).^2));
    vs=spline_matrix(mup,mu);
    v=[vs,zeros(n,n-np)];
else
    % to a denser media
    critical=cos(asin(sqrt(epsilon_i/epsilon_t)));  
    np=length(find(mu>critical));
    mup=zeros(1,n);
    mup=sqrt(1-epsilon_i/epsilon_t*(1-mu.^2));
     vs=spline_matrix(mup,mu);
     v=[vs(1:np,:);zeros(n-np,n)];
 end
 