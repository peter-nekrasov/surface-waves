function v = get_rhs_vec_flex_pt_src(k,nu,n,R,src)

if nargin < 5
    src = [0; 0];
end

rho = sqrt(src(1)^2 + src(2)^2);
theta = atan2(src(2),src(1));

hn = besselj(n,k*rho)*besselh(n,k*R);
hnp = besselj(n,k*rho)*k/2*(besselh(n-1,k*R)-besselh(n+1,k*R));
hnpp = besselj(n,k*rho)*(-k*R*besselh(n-1,k*R)...
        +(n+n^2-k^2*R^2)*besselh(n,k*R))/R^2;
hnppp = besselj(n,k*rho)*(k*R*(2+n^2-k^2*R^2)*besselh(-1+n,k*R)...
         - (1+n)*(2*n+n^2-k^2*R^2)*besselh(n,k*R))/R^3;

kn = besseli(n,k*rho)*besselk(n,k*R);
knp = besseli(n,k*rho)*k/2*(-besselk(n-1,k*R)-besselk(n+1,k*R));
knpp = besseli(n,k*rho)*(k*R*besselk(n-1,k*R)...
        +(n+n^2+k^2*R^2)*besselk(n,k*R))/R^2;
knppp = besseli(n,k*rho)*(-k*R*(2+n^2+k^2*R^2)*besselk(-1+n,k*R)...
         - (1+n)*(2*n+n^2+k^2*R^2)*besselk(n,k*R))/R^3;

gn = hn*1i/(8*k^2) - kn/(4*pi*k^2);
gnp = hnp*1i/(8*k^2) - knp/(4*pi*k^2);
gnpp = hnpp*1i/(8*k^2) - knpp/(4*pi*k^2);
gnppp = hnppp*1i/(8*k^2) - knppp/(4*pi*k^2);


v = [gnpp + nu/R*gnp - nu/R^2*n^2*gn; 
    gnppp + 1/R*gnpp - (1+n^2*(2-nu))/R^2*gnp + (3-nu)/R^3*n^2*gn];

v(1) = v(1)*exp(-1i*n*theta);
v(2) = v(2)*exp(-1i*n*theta);

end