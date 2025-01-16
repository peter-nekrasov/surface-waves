clear 
close all
addpath(genpath('../..'))

R = sqrt(2);
nu = 0.3;
n = 0;
gamma = 0.01;
beta = 5+1i;
k = beta^(1/4);
w = 1;

f0 =  get_rhs_vec_flex_pt_src(k,nu,n,R);
A = get_lhs_for_bcs(k,nu,n,R);
epsilon = 1E-8;
d = 1;

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.ta = 0;
cparams.tb = R;

fcurve = @(t) [t(:).';0*t(:).'];
int_chnkr = chunkerfunc(fcurve,cparams);
int_chnkr = sort(int_chnkr);

cparams.maxchunklen = 1 / abs(k);
cparams.ta = R;
cparams.tb = 100; % 10E5
ext_chnkr = chunkerfunc(fcurve,cparams);
ext_chnkr = sort(ext_chnkr);

chnkr = merge([int_chnkr ext_chnkr]);
chnkr = sort(chnkr);

figure(1)
clf
plot(int_chnkr,'.')
hold on
plot(ext_chnkr,'.')

f_int = flex2d.hkdiffgreen(k,[0;0],int_chnkr.r) / (2*k^2);
f_ext = flex2d.hkdiffgreen(k,[0;0],ext_chnkr.r) / (2*k^2);
f_tot = flex2d.hkdiffgreen(k,[0;0],chnkr.r) / (2*k^2);

sig_int_n = zeros(int_chnkr.npt,1);
sig_ext_n = zeros(ext_chnkr.npt,1);
sig_tot_n = [sig_int_n; sig_ext_n];

% build the desired kernels

gkern =  @(s,t) flex2d.modal_gf(k, s, t, 0); 
skern =  @(s,t) axissymlap2d.kern(s, t, [0;0], 's') / (2*pi^2); 

opts = [];
opts.sing = 'log';

G = chunkermat(chnkr,gkern, opts);
S = chunkermat(chnkr,skern, opts);
GS = G*S;

figure(1)
title('\sigma')

while d > epsilon

    plot(chnkr.r(1,:),real(sig_tot_n))
    hold on

    GSu = -gamma*GS*(sig_tot_n+f_tot);
    
    sig_ext_np12 = GSu(int_chnkr.npt+1:end);
    sig_tot_np12 = [sig_int_n; sig_ext_np12];

    Su = S*(sig_tot_np12+f_tot);
    sig_int_np1 = - f_int + w^2*Su(1:int_chnkr.npt)/2;

    dmat = lege.dermat(chnkr.k);
    dd = ext_chnkr.d(1,:,1);

    gsup = dmat*GSu(int_chnkr.npt+1:int_chnkr.npt+16)./dd(:);
    gsupp = dmat*gsup./dd(:);
    gsuppp = dmat*gsupp./dd(:);

    b = - [gsupp(1) + nu/R*gsup(1); gsuppp(1) + 1/R*gsupp(1) - 1/R^2*gsup(1)] - f0;
    v = A \ b;

    sig_ext_np1 = sig_ext_np12 + v(1)*besselh(n,k*ext_chnkr.r(1,:)).' + ...
        v(2)*besselk(n,k*ext_chnkr.r(1,:)).';

    sig_tot_np1 = [sig_int_np1; sig_ext_np1];

    d = max(abs(sig_tot_n - sig_tot_np1)) / max(abs(sig_tot_np1))

    sig_tot_n = sig_tot_np1;
    sig_int_n = sig_int_np1;
    sig_ext_n = sig_ext_np1;

end

phi_n_tot = sig_tot_n + f_tot;
phi_tot = S*phi_n_tot;

figure(2)
tiledlayout(1,2);
title('\phi_n')
nexttile
plot(chnkr.r(1,:),real(phi_n_tot))
nexttile
plot(chnkr.r(1,:),real(phi_tot))

return 

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


