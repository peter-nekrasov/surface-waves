clear 
close all
addpath(genpath('../..'))

R = sqrt(2);
k = 5+0.2i;
nu = 0.3;
n = 0;


f0 =  get_rhs_vec_flex_pt_src(k,nu,n,R);
A = get_lhs_for_bcs(k,nu,n,R);
epsilon = 1E-9;
d = 1;

% create domains in chunkie 
cparams = [];
cparams.maxchunklen = 8 / abs(k);
cparams.ifclosed = false;
cparams.ta = 0;
cparams.tb = R;

fcurve = @(t) [t(:).';0*t(:).'];
int_chnkr = chunkerfunc(fcurve,cparams);

cparams.ta = R;
cparams.tb = 100; % 10E5
ext_chnkr = chunkerfunc(fcurve,cparams);

chnkr = merge([int_chnkr ext_chnkr]);

figure(1)
clf
plot(int_chnkr,'.')
hold on
plot(ext_chnkr,'.')

sig_int = zeros(int_chnkr.npt,1);
sig_ext = zeros(ext_chnkr.npt,1);

f_int = hkdiffgreen(k,[0;0],int_chnkr.r);
f_ext = hkdiffgreen(k,[0;0],ext_chnkr.r);
f0 = hkdiffgreen(k,[0;0],chnkr.r);


fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'modal-gf', coefs);        % build the desired kernel


return 

while d > epsilon
    j = 0;
end

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


