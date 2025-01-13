clear 
close all
addpath(genpath('../..'))

R = 1;
k = 5;
nu = 0.3;
n = 0;

v =  get_rhs_vec_flex_pt_src(k,nu,n,R);
A = get_lhs_for_bcs(k,nu,n,R);

coefs = A \ v;

src = [0;0];
targ = [(0.1:0.1:10);(0.1:0.1:10)*0];

[h0,~] = helmdiffgreen(k,src,targ);
[k0,~] = helmdiffgreen(1i*k,src,targ);
[g0,~] = hkdiffgreen(k,src,targ);

sol = coefs(1)*h0 + coefs(2)*k0;

err = max(sol - g0)

rmpath(genpath('../..'))


