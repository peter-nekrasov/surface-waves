clear 
close all
addpath(genpath('../..'))

R = 3*sqrt(5);
nu = 0.3;
n = 0;
k = 2+0.02i;
s = 2;

w = 1;
g = 1;
alpha = 1;
gamma = 1;
beta = 1;

% f0 = get_rhs_vec_gaussian(s,nu,n,R);
% A = get_lhs_for_bcs(k,nu,n,R);
% vf = A \ f0; 
% vf = [0;0];
% epsilon = 1E-8;
% d = 1;

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = R;

fcurve = @(t) [t(:).';0*t(:).'];

cparams.ta = 0;
cparams.tb = R; % 10E5

pref = [];
pref.nchmax = 2;
int_chnkr = chunkerfunc(fcurve,cparams,pref);
int_chnkr = sort(int_chnkr);

pref = [];
pref.nchmax = 199/R;
cparams.maxchunklen = 1.1*R;
cparams.ta = R;
cparams.tb = 300; 

ext_chnkr = chunkerfunc(fcurve,cparams);
ext_chnkr = sort(ext_chnkr);

chnkr = merge([int_chnkr ext_chnkr]);

 
figure(1)
plot(int_chnkr,'.')
hold on 
plot(ext_chnkr,'.')
title('chnkrs')


f_ext = exp(-(ext_chnkr.r(1,:)-8*R).^2/(2*s^2)); % + vf(1)*besselh(n,k*chnkr.r(1,:)) + vf(2)*besselk(n,k*chnkr.r(1,:));
f_ext = f_ext.';

f_int = 0*int_chnkr.r(1,:).';

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

figure(2);
plot(int_chnkr.r(1,:),real(f_int),ext_chnkr.r(1,:),real(f_ext))
title('f')

sig_n = zeros(chnkr.npt,1);

% build the desired kernels

gkern =  @(s,t) flex2d.ext_modal_gf(k,R,nu,s,t,0); 
skern =  @(s,t) axissymlap2d.kern(s, t, [0;0], 's') / (2*pi^2); 

opts = [];
opts.sing = 'log';

G = chunkermat(ext_chnkr,gkern, opts);
S = chunkermat(chnkr,skern, opts);

% make derivatives 

lhs = eye(chnkr.npt);
lhs(1:int_chnkr.npt,:) = lhs(1:int_chnkr.npt,:) - w^2/g*S(1:int_chnkr.npt,:);
lhs(int_chnkr.npt+1:end,:) = lhs(int_chnkr.npt+1:end,:) + gamma/alpha*G*S(int_chnkr.npt+1:end,:);
rhs = [f_int; G*f_ext];
sol = lhs\rhs;

figure(3)
plot(chnkr.r(1,:),real(sol),chnkr.r(1,:),imag(sol))
title('solution ($\hat{\sigma}$)','Interpreter','latex')
legend('real','imaginary')
xlim([0 100])

Su = S*sol(:);
uint = sol(1:int_chnkr.npt);
Suint = Su(1:int_chnkr.npt);
Suext = Su(int_chnkr.npt+1:end);
dmat = lege.dermat(16);
sol = squeeze(reshape(sol,size(chnkr.r(1,:,:))));
dd = squeeze(ext_chnkr.d(1,:,:));
u = sol(:,3:end);
dudr = dmat*u./dd;
d2udr2 = dmat*dudr./dd;
d3udr3 = dmat*d2udr2./dd;
d4udr4 = dmat*d3udr3./dd;

%hold on
%plot(chnkr.r(1,:),real(dudr(:)),chnkr.r(1,:),imag(dudr(:)))

rs = ext_chnkr.r(1,:);
f_recon = d4udr4(:)+2./rs(:).*d3udr3(:)-1./rs(:).^2.*d2udr2(:)+1./rs(:).^3.*dudr(:)+Suext(:)-k^4*u(:);

figure(4)
plot(rs,real(f_recon-f_ext),rs,imag(f_recon-f_ext))
title('Lu + S[u](= f)')

figure(5)
rs = int_chnkr.r(1,:);
plot(rs,real(uint-w^2/g*Suint),rs,imag(uint-w^2/g*Suint))
title('u + S[u](= 0)')

% check that BCs are satisfied by the solution

u1 = u(:,1);
ds = ext_chnkr.d(1,:,1).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);

ddr = pols.'*v2c*dmat*(u1./ds);
d2dr2 = pols.'*v2c*dmat^2*(u1./ds.^2);
d3dr3 = pols.'*v2c*dmat^3*(u1./ds.^3);

bc1 = d2dr2 + nu/R*ddr; 
bc2 = d3dr3 + 1/R*d2dr2 - 1/R^2*ddr;

rmpath(genpath('../..'))
