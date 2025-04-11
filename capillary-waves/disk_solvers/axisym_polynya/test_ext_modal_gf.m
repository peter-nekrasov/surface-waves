%% Plotting modal green's function vs exterior modal green's function

clear 
close all
addpath(genpath('../..'))

R = sqrt(2);
rs = R:R/10:60*R;
k = 2; %+0.05i;

src = [];
src.r = [8; 0];

targ = [];
targ.r = [rs; 0*rs];

nu = 0.3;

g_ext = flex2d.ext_modal_gf(k,R,nu, src, targ, 0);
g = flex2d.modal_gf(k, src, targ, 0);

% compare bounded modal GF to regular modal GF

figure(1);
plot(rs,real(g_ext),rs,real(g))
legend('exterior modal gf', 'modal gf')


%% Checking that the modal green's function satisfies the radial ODE

[rs,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);
rs = rs/4 + 0.25 + R;
dd = ones(16,1)/4;

src = [];
src.r = [5; 0];

targ = [];
targ.r = [rs.'; 0*rs.'];
nu = 0.3;

g_ext = flex2d.ext_modal_gf(k,R,nu, src, targ, 0);
g = flex2d.modal_gf(k, src, targ, 0);

% compare bounded modal GF to regular modal GF

figure(1);
plot(rs,real(g_ext),rs,real(g))
legend('ext modal gf', 'modal gf')

% g = flex2d.hkdiffgreen(k,src.r,targ.r);
dmat = lege.dermat(16);
dgdr = dmat*g_ext./dd;
d2gdr2 = dmat*dgdr./dd;
d3gdr3 = dmat*d2gdr2./dd;
d4gdr4 = dmat*d3gdr3./dd;

err = d4gdr4+2./rs.*d3gdr3-1./rs.^2.*d2gdr2+1./rs.^3.*dgdr-k^4*g_ext;

% check that BCs are satisfied 

ddr = pols.'*v2c*dmat*(g_ext./dd);
d2dr2 = pols.'*v2c*dmat^2*(g_ext./dd.^2);
d3dr3 = pols.'*v2c*dmat^3*(g_ext./dd.^3);

bc1 = d2dr2 + nu/R*ddr; 
bc2 = d3dr3 + 1/R*d2dr2 - 1/R^2*ddr;

%% Checking that \int G f gives the correct solution to the ODE

k = 2+0.5i;
R = sqrt(2);

cparams = [];
cparams.ifclosed = false;

fcurve = @(t) [t(:).';0*t(:).'];
cparams.maxchunklen = 4 / abs(k);
cparams.ta = R;
cparams.tb = 20; % 10E5
chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

f = (chnkr.r(1,:)-5).'.*normpdf(5*chnkr.r(1,:)-25).';

figure(1);
t = tiledlayout(1,3);
nexttile 
plot(chnkr.r(1,:),f(:))
title('f(r)')

gkern =  @(s,t) flex2d.ext_modal_gf(k,R,nu, s, t, 0);

opts = [];
opts.sing = 'log';
G = chunkermat(chnkr,gkern, opts);

sol = G*f;
nexttile 
plot(chnkr.r(1,:),real(sol),chnkr.r(1,:),imag(sol))
title('u = \int G f')
legend('real part','imaginary part')
xlim([0 10])


dmat = lege.dermat(16);
sol = squeeze(reshape(sol,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));

u = sol;
dudr = dmat*u./dd;
d2udr2 = dmat*dudr./dd;
d3udr3 = dmat*d2udr2./dd;
d4udr4 = dmat*d3udr3./dd;

%hold on
%plot(chnkr.r(1,:),real(dudr(:)),chnkr.r(1,:),imag(dudr(:)))

rs = chnkr.r(1,:);
f_recon = d4udr4(:)+2./rs(:).*d3udr3(:)-1./rs(:).^2.*d2udr2(:)+1./rs(:).^3.*dudr(:)-k^4*u(:);

nexttile
plot(rs,real(f_recon),rs,imag(f_recon))
title('L u (= f)')
xlim([0 20])


% check that BCs are satisfied by the solution

u1 = u(:,1);
ds = chnkr.d(1,:,1).';

ddr = pols.'*v2c*dmat*(u1./ds);
d2dr2 = pols.'*v2c*dmat^2*(u1./ds.^2);
d3dr3 = pols.'*v2c*dmat^3*(u1./ds.^3);

bc1 = d2dr2 + nu/R*ddr; 
bc2 = d3dr3 + 1/R*d2dr2 - 1/R^2*ddr;



