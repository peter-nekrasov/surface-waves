clear 
close all
addpath(genpath('../..'))

R = sqrt(5);
nu = 0.3;
n = 0;
gamma = 1;
k = 2+0.1i;
s = 2;

% f0 = get_rhs_vec_gaussian(s,nu,n,R);
% A = get_lhs_for_bcs(k,nu,n,R);
% vf = A \ f0; 
% vf = [0;0];
epsilon = 1E-8;
d = 1;

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 8 / abs(k);
cparams.ta = R;
cparams.tb = 300; % 10E5

fcurve = @(t) [t(:).';0*t(:).'];

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

% figure(1)
% clf
% plot(chnkr,'.')

f = exp(-(chnkr.r(1,:)-10*R).^2/(2*s^2)); % + vf(1)*besselh(n,k*chnkr.r(1,:)) + vf(2)*besselk(n,k*chnkr.r(1,:));
f = f.';

% [~,~,u] = lege.exps(16);
% [pols,~] = lege.pols(-1,15);

% dmat = lege.dermat(chnkr.k);
% dd = chnkr.d(1,:,1);
% 
% gsup = dmat*f(1:16)./dd(:);
% gsupp = dmat*gsup./dd(:);
% gsuppp = dmat*gsupp./dd(:);
% 
% gsup = pols.'*u*gsup;
% gsupp = pols.'*u*gsupp;
% gsuppp = pols.'*u*gsuppp;
% 
% b = [gsupp + nu/R*gsup; gsuppp + 1/R*gsupp - 1/R^2*gsup];

figure(1);
plot(chnkr.r(1,:),real(f),chnkr.r(1,:),imag(f))
title('f')
hold on

sig_n = zeros(chnkr.npt,1);

% build the desired kernels

gkern =  @(s,t) flex2d.ext_modal_gf(k,R,nu, s, t, 0); 
skern =  @(s,t) axissymlap2d.kern(s, t, [0;0], 's') / (2*pi^2); 

opts = [];
opts.sing = 'log';

G = chunkermat(chnkr,gkern, opts);
S = chunkermat(chnkr,skern, opts);

% make derivatives 

lhs = eye(chnkr.npt) + gamma*G*S;
rhs = G*f;
sol = lhs\rhs;

figure(2)
plot(chnkr.r(1,:),real(sol),chnkr.r(1,:),imag(sol))
title('solution')

dmat = lege.dermat(16);
sol = squeeze(reshape(sol,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
u = sol;
dudr = dmat*u./dd;
d2udr2 = dmat*dudr./dd;
d3udr3 = dmat*d2udr2./dd;
d4udr4 = dmat*d3udr3./dd;
Su = S*sol(:);

%hold on
%plot(chnkr.r(1,:),real(dudr(:)),chnkr.r(1,:),imag(dudr(:)))

rs = chnkr.r(1,:);
f_recon = d4udr4(:)+2./rs(:).*d3udr3(:)-1./rs(:).^2.*d2udr2(:)+1./rs(:).^3.*dudr(:)+Su(:)-k^4*u(:);

figure(3)
plot(rs,real(f_recon),rs,imag(f_recon))
title('Lu + S[u](= f)')

figure(4)
plot(rs,real(f_recon-f),rs,imag(f_recon-f))
title('residual error')

% check that BCs are satisfied by the solution

u1 = u(:,1);
ds = chnkr.d(1,:,1).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);

ddr = pols.'*v2c*dmat*(u1./ds);
d2dr2 = pols.'*v2c*dmat^2*(u1./ds.^2);
d3dr3 = pols.'*v2c*dmat^3*(u1./ds.^3);

bc1 = d2dr2 + nu/R*ddr; 
bc2 = d3dr3 + 1/R*d2dr2 - 1/R^2*ddr;

rmpath(genpath('../..'))


