function A = get_lhs_for_bcs(k,nu,n,R)


hn = besselh(n,k*R);
hnp = k/2*(besselh(n-1,k*R)-besselh(n+1,k*R));
hnpp = (-k*R*besselh(n-1,k*R)...
        +(n+n^2-k^2*R^2)*besselh(n,k*R))/R^2;
hnppp = (k*R*(2+n^2-k^2*R^2)*besselh(-1+n,k*R)...
         - (1+n)*(2*n+n^2-k^2*R^2)*besselh(n,k*R))/R^3;

kn = besselk(n,k*R);
knp = k/2*(-besselk(n-1,k*R)-besselk(n+1,k*R));
knpp = (k*R*besselk(n-1,k*R)...
        +(n+n^2+k^2*R^2)*besselk(n,k*R))/R^2;
knppp = (-k*R*(2+n^2+k^2*R^2)*besselk(-1+n,k*R)...
         - (1+n)*(2*n+n^2+k^2*R^2)*besselk(n,k*R))/R^3;

A = [hnpp + nu/R*hnp - nu/R^2*n^2*hn, knpp + nu/R*knp - nu/R^2*n^2*kn; 
    hnppp + 1/R*hnpp - (1+n^2*(2-nu))/R^2*hnp + (3-nu)/R^3*n^2*hn, ...
    knppp + 1/R*knpp - (1+n^2*(2-nu))/R^2*knp + (3-nu)/R^3*n^2*kn];

end