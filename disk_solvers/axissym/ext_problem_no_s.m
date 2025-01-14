clear 
close all
addpath(genpath('../..'))

R = 1.6;
k = 5;
nu = 0.3;
n = 0;

v =  get_rhs_vec_flex_pt_src(k,nu,n,R);
A = get_lhs_for_bcs(k,nu,n,R);

coefs = A \ v;

src = [0;0];
targ = [(R:R/100:5*R);(R:R/100:5*R)*0];

dr = sqrt((targ(1,:)-src(1)).^2 + (targ(2,:)-src(2)).^2);

h0 = besselh(0,k*dr);
k0 = besselk(0,k*dr);
g0 = h0*1i/(8*k^2) - k0/(4*pi*k^2);

sol = coefs(1)*h0 + coefs(2)*k0;

err = max(abs(sol - g0)) / max(abs(g0))

rmpath(genpath('../..'))


