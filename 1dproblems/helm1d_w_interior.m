clear 
close all
addpath(genpath('..'))

L = 5;
k = 1+0.2i;

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 2 / abs(k);
cparams.ta = -100;
cparams.tb = -L; 

fcurve = @(t) [t(:).'; 0*t(:).'];

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

cparams.ta = -L;
cparams.tb = L; 

int_chnkr = chunkerfunc(fcurve,cparams);
int_chnkr = sort(int_chnkr);

tot_chnkr = merge([chnkr int_chnkr]);
tot_chnkr = sort(tot_chnkr);

figure(1)
plot(chnkr,'.')
hold on
plot(int_chnkr,'.')
title('chnkr')

% s = 2;
% f = -(chnkr.r(1,:)+8*L).*exp(-(chnkr.r(1,:)+8*L).^2/(2*s^2)); % + vf(1)*besselh(n,k*chnkr.r(1,:)) + vf(2)*besselk(n,k*chnkr.r(1,:));
% f = f.';

% figure(2)
% plot(chnkr.r(1,:), f)
% title('f (RHS)')

% method 1, solve exterior problem, then interior integral equation

bdry = [];
bdry.r = [-L; 0];

L0 = -10*L; 
src1 = []; src1.r = [L0; 0]; 

[~,grad] = helm1d.green(k,src1,bdry);

mu1 = 2*grad;

gkern =  @(s,t) helm1d.green(k,s,t); 
skern =  @(s,t) lap2d.kern(s,t,'s'); 

opts = [];
opts.sing = 'log';

% G = chunkermat(chnkr,gkern, opts);
S = chunkermat(tot_chnkr,skern, opts);

uext = gkern(src1,chnkr) + mu1*gkern(bdry,chnkr);

N = int_chnkr.npt;

rhs = S(end-N+1:end,1:end-N)*uext;
lhs = eye(N) - S(end-N+1:end,end-N+1:end);

uint = lhs \ rhs;
usol = [uext; uint];

figure(3)
plot(tot_chnkr.r(1,:),real(usol),tot_chnkr.r(1,:),imag(usol))
title('u (solution)')


% checking the solution

dmat = lege.dermat(16);
u = squeeze(reshape(uext,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = chnkr.r(1,:);
f_recon = d2udx2(:)+k^2*u(:);

figure(4)
plot(rs,real(f_recon),rs,imag(f_recon))
title('Exterior error (Lu - f)')

Su = S*usol;
int_err = usol(end-N+1:end) - Su(end-N+1:end);
figure(5)
plot(int_chnkr.r(1,:), real(int_err), int_chnkr.r(1,:), imag(int_err))
title('Interior error (u - S[u])')

% check that the BC is satisfied by the solution

u1 = u(:,end);
ds = chnkr.d(1,:,end).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(1,15);
dudn = pols.'*v2c*dmat*(u1./ds);


% method 2, solve first kind integral equation using adjoint density

S = chunkermat(int_chnkr,skern, opts);
G = chunkermat(int_chnkr,gkern, opts);

lhs = G - S*G;
%lhs(1:N,N+1) = gkern(bdry,int_chnkr) - S*gkern(bdry,int_chnkr);

[~,grad] = helm1d.green(k,int_chnkr,bdry);
%lhs(N+1,1:N) =  -grad.*int_chnkr.wts(:).';
%lhs(N+1,N+1) = 1/2;

S = chunkermat(tot_chnkr,skern, opts);
rhs = - gkern(src1,int_chnkr) + S(end-N+1:end,:)*gkern(src1,tot_chnkr); %; grad2];

mu = lhs \ rhs;

mu1 = mu(end);
mu = mu(1:N);

G = chunkermat(tot_chnkr,gkern, opts);
usol = G(:,end-N+1:end)*mu + gkern(src1,tot_chnkr); % + mu1*gkern(bdry,tot_chnkr);

figure(6)
plot(int_chnkr.r(1,:),real(mu),int_chnkr.r(1,:),imag(mu))
title('\mu (interior density)')


figure(7)
plot(tot_chnkr.r(1,:),real(usol),tot_chnkr.r(1,:),imag(usol))
title('u (solution)')

return

% method 3, solve preconditioned integral equation using adjoint density

lhs = zeros(N+1);
S = chunkermat(int_chnkr,skern, opts);

lhs(1:N,1:N) = eye(N) - S;
lhs(1:N,N+1) = -skern(bdry,int_chnkr);

[~,grad] = helm1d.green(k,int_chnkr,bdry);

lhs(N+1,1:N) = -(grad(:).*int_chnkr.wts(:)).';
lhs(N+1,N+1) = 1/2;

[~,grad] = helm1d.green(k,src1,bdry);

rhs = [skern(src1,int_chnkr); grad];

mu = lhs \ rhs;
mu1 = mu(end);
mu = mu(1:N);


G = chunkermat(tot_chnkr,gkern, opts);
usol = G(:,end-N+1:end)*mu + gkern(src1,tot_chnkr) + mu1*gkern(bdry,tot_chnkr);

figure(6)
plot(int_chnkr.r(1,:),real(mu),int_chnkr.r(1,:),imag(mu))
title('\mu (interior density)')


figure(7)
plot(tot_chnkr.r(1,:),real(usol),tot_chnkr.r(1,:),imag(usol))
title('u (solution)')

return

% checking the solution

dmat = lege.dermat(16);
u = squeeze(reshape(uext,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = chnkr.r(1,:);
f_recon = d2udx2(:)+k^2*u(:);

figure(4)
plot(rs,real(f_recon),rs,imag(f_recon))
title('Lu (= f)')

err = (f_recon(:) - f(:))/max(f(:)).*wts(:).^2;

figure(5)
plot(rs,real(err),rs,imag(err))
title('absolute error (residual)')

% check that the BC is satisfied by the solution

u1 = u(:,end);
ds = chnkr.d(1,:,end).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(1,15);
dudn = pols.'*v2c*dmat*(u1./ds);

rmpath(genpath('..'))