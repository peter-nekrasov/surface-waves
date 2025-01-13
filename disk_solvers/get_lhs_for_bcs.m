function A = get_lhs_for_bcs(k,nu,n,R)

src = [0; 0];
targ = [R; 0];

[h0,h0p,h0pp,h0ppp] = helmdiffgreen(k,src,targ);
[k0,k0p,k0pp,k0ppp] = helmdiffgreen(1i*k,src,targ);

h0p = h0p(1,1,1);
h0pp = h0pp(1,1,1);
h0ppp = h0ppp(1,1,1);

k0p = k0p(1,1,1);
k0pp = k0pp(1,1,1);
k0ppp = k0ppp(1,1,1);

A = [h0pp + nu/R*h0p - nu/R^2*n^2*h0, k0pp + nu/R*k0p - nu/R^2*n^2*k0; 
    h0ppp + 1/R*h0pp - (1+n^2*(2-nu))/R^2*h0p + (3-nu)/R^3*n^2*h0, ...
    k0ppp + 1/R*k0pp - (1+n^2*(2-nu))/R^2*k0p + (3-nu)/R^3*n^2*k0];

end