clear 
close all
addpath(genpath('../..'))

% Parameters
R = 1.3; % size of disk
k = 5; % wavenumber
nu = 0.3; % poisson ratio
N = 20; % number of modes 
src_loc = [0.2; 0.1]; % location of pt src

% get true solution 
targ = [(R:R/50:10*R);(R:R/50:10*R)*0];
dr = sqrt((targ(1,:)-src_loc(1)).^2 + (targ(2,:)-src_loc(2)).^2);
h0 = besselh(0,k*dr);
k0 = besselk(0,k*dr);
g0 = h0*1i/(8*k^2) - k0/(4*pi*k^2);

sol = 0;

dr = sqrt((targ(1,:)).^2 + (targ(2,:)).^2);

for n = -N:N
    v = get_rhs_vec_flex_pt_src(k,nu,n,R,src_loc);
    A = get_lhs_for_bcs(k,nu,n,R);
    coefs = A \ v;

    hn = besselh(n,k*dr);
    kn = besselk(n,k*dr);

    sol = sol + coefs(1)*hn + coefs(2)*kn;
end

plot(dr,real(sol))

err = max(abs(sol - g0)) / max(abs(g0))

rmpath(genpath('../..'))






