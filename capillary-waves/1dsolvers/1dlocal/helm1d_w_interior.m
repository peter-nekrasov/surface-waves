L = 10;
k = 1+0.05i;

left_bd = 0;
right_bd = L;

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 2 / abs(k);
cparams.ta = -200;
cparams.tb = left_bd; 

fcurve = @(t) [t(:).'; 0*t(:).'];

left_chnkr = chunkerfunc(fcurve,cparams);
left_chnkr = sort(left_chnkr);

cparams.ta = right_bd;
cparams.tb = 200; 

right_chnkr = chunkerfunc(fcurve,cparams);
right_chnkr = sort(right_chnkr);

cparams.ta = left_bd;
cparams.tb = right_bd; 
% cparams.chsmall = L/50;

cparams.maxchunklen = 2 / abs(k);
int_chnkr = chunkerfunc(fcurve,cparams);
int_chnkr = sort(int_chnkr);

tot_chnkr = merge([left_chnkr int_chnkr right_chnkr]);
tot_chnkr = sort(tot_chnkr);

N_int = int_chnkr.npt;
N_left = left_chnkr.npt;
N_right = right_chnkr.npt;

figure(1)
plot(left_chnkr,'.')
hold on
plot(int_chnkr,'.')
hold on 
plot(right_chnkr,'.')
title('chnkr')

L0 = -5*L; 

s = 2;
f = exp(-(left_chnkr.r(1,:)-L0).^2/(2*s^2));
f = f.';

figure(2)
plot(left_chnkr.r(1,:), f)
title('f (RHS)')

%% method 1, solve exterior problem, then solve interior equation

bdry1 = [];
bdry1.r = [left_bd; 0];

bdry2 = [];
bdry2.r = [right_bd; 0];

src1 = []; src1.r = [L0; 0]; 

[val,grad] = helm1d.green(k,left_chnkr,bdry1);
grad = grad(:,:,1);

wts = left_chnkr.wts;

mu1 = 2*grad*(f.*wts(:));

gkern =  @(s,t) helm1d.green(k,s,t); 
skern =  @(s,t) chnk.lap2d.kern(s,t,'s'); 

opts = [];
opts.sing = 'removable';

opts2 = [];
opts2.sing = 'log';

G = chunkermat(left_chnkr,gkern, opts);
S = chunkermat(tot_chnkr,skern, opts2);

uleft = G*f + mu1*gkern(bdry1,left_chnkr);
uright = 0*right_chnkr.r(1,:).';

rhs = S(N_left+1:N_left+N_int,1:N_left)*uleft;
lhs = eye(N_int) - S(N_left+1:N_left+N_int,N_left+1:N_left+N_int);

uint = lhs \ rhs;
usol = [uleft; uint; uright];

usol_analytic = usol;

figure(3)
plot(tot_chnkr.r(1,:),real(usol),tot_chnkr.r(1,:),imag(usol))
title('u (solution)')

% checking the solution

dmat = lege.dermat(16);
u = squeeze(reshape(uleft,size(left_chnkr.r(1,:,:))));
dd = squeeze(left_chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = left_chnkr.r(1,:);
f_recon = d2udx2(:)+k^2*u(:);

figure(4)
plot(rs,real(f-f_recon),rs,imag(f-f_recon))
title('Exterior error (Lu - f)')

Su = S*usol;
int_err = usol(N_left+1:N_left+N_int) - Su(N_left+1:N_left+N_int);
figure(5)
plot(int_chnkr.r(1,:), real(int_err), int_chnkr.r(1,:), imag(int_err))
title('Interior error (u - S[u])')

% check that the BC is satisfied by the solution

u1 = u(:,end);
ds = left_chnkr.d(1,:,end).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(1,15);
dudn = pols.'*v2c*dmat*(u1./ds);

%% method 2, solve first kind integral equation using adjoint formulation

lhs = zeros(N_int+2,N_int+2);

% dd = squeeze(int_chnkr.d(1,:,1));
% lhs(N_int+3,1:16) = dmat(1,:)./dd;
% dd = squeeze(int_chnkr.d(1,:,end));
% lhs(N_int+4,end-17:end-2) =  dmat(end,:)./dd;

tot_chnkr = merge([left_chnkr int_chnkr right_chnkr]);
tot_chnkr = sort(tot_chnkr);

S = chunkermat(tot_chnkr,skern, opts2);
G = chunkermat(tot_chnkr,gkern, opts);

lhs(1:N_int,1:N_int) = chunkermat(int_chnkr,gkern, opts) - S(N_left+1:N_left+N_int,:)*G(:,N_left+1:N_left+N_int);

lhs(1:N_int,N_int+1) = gkern(bdry1,int_chnkr) - S(N_left+1:N_left+N_int,:)*gkern(bdry1,tot_chnkr);
lhs(1:N_int,N_int+2) = gkern(bdry2,int_chnkr) - S(N_left+1:N_left+N_int,:)*gkern(bdry2,tot_chnkr);

[~,grad] = helm1d.green(k,int_chnkr,bdry1);
grad = grad(:,:,1);
lhs(N_int+1,1:N_int) = grad.*int_chnkr.wts(:).';
lhs(N_int+1,N_int+1) = -1/2;

[~,grad] = helm1d.green(k,bdry2,bdry1);
grad = grad(:,:,1);
lhs(N_int+1,N_int+2) = grad;

[~,grad] = helm1d.green(k,int_chnkr,bdry2);
grad = grad(:,:,1);
lhs(N_int+2,1:N_int) = grad.*int_chnkr.wts(:).';
lhs(N_int+2,N_int+2) = 1/2;

[~,grad] = helm1d.green(k,bdry1,bdry2);
grad = grad(:,:,1);
lhs(N_int+2,N_int+1) = grad;

[~,grad2] = helm1d.green(k,left_chnkr,bdry1);
grad2 = grad2(:,:,1);

[~,grad3] = helm1d.green(k,left_chnkr,bdry2);
grad3 = grad3(:,:,1);

wts = left_chnkr.wts;

rhs = [- G(N_left+1:N_left+N_int,1:N_left)*f + S(N_left+1:N_left+N_int,:)*G(:,1:N_left)*f ; -grad2*(f.*wts(:)); -grad3*(f.*wts(:))]; % 0; 0];

mu = lhs \ rhs;
mu2 = mu(end);
mu1 = mu(end-1);
mu = mu(1:N_int);

mu_true = mu;
mu1_true = mu1;
mu2_true = mu2;

G = chunkermat(tot_chnkr,gkern, opts);
usol = G(:,N_left+1:N_left+N_int)*mu + G(:,1:N_left)*f + mu1*gkern(bdry1,tot_chnkr) + mu2*gkern(bdry2,tot_chnkr);

usol_fk = usol;

figure(6)
plot(int_chnkr.r(1,:),real(mu),int_chnkr.r(1,:),imag(mu))
title('\mu (interior density)')

figure(7)
plot(tot_chnkr.r(1,:),real(usol),tot_chnkr.r(1,:),imag(usol))
title('u (solution)')

% figure(8)
% loglog(int_chnkr.r(1,:),abs(mu))
% hold on
% loglog(int_chnkr.r(1,:),10*int_chnkr.r(1,:).^(-1.8))
% legend('mu','r^{-1.8}')


% figure(8)
% plot(tot_chnkr.r(1,:),real(usol - usol_analytic),tot_chnkr.r(1,:),imag(usol - usol_analytic))
% title('u (solution)')

% checking the solution

dmat = lege.dermat(16);
uext = usol(1:N_left);
u = squeeze(reshape(uext,size(left_chnkr.r(1,:,:))));
dd = squeeze(left_chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = left_chnkr.r(1,:);
f_recon = d2udx2(:)+k^2*u(:);

figure(8)
plot(rs,real(f_recon-f),rs,imag(f_recon-f))
title('Exterior error (Lu - f)')

Su = S*usol;
int_err = usol(N_left+1:N_left+N_int) - Su(N_left+1:N_left+N_int);

figure(9)
plot(int_chnkr.r(1,:), real(int_err), int_chnkr.r(1,:), imag(int_err))
title('Interior error (u - S[u])')

% check that the BC is satisfied by the solution

u1 = u(:,end);
ds = left_chnkr.d(1,:,end).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(1,15);
dudn = pols.'*v2c*dmat*(u1./ds);

%% method 3, solve preconditioned integral equation . doesn't work 

%{
lhs = zeros(N_int+2);

S = chunkermat(int_chnkr,skern, opts2);

lhs(1:N_int,1:N_int) = eye(N_int) - S;

lhs(1:N_int,N_int+1) = -skern(bdry1,int_chnkr);
lhs(1:N_int,N_int+2) = -skern(bdry2,int_chnkr);

[~,grad] = helm1d.green(k,int_chnkr,bdry1);
grad = grad(:,:,1);
lhs(N_int+1,1:N_int) = grad.*int_chnkr.wts(:).';
lhs(N_int+1,N_int+1) = -1/2;

[~,grad] = helm1d.green(k,bdry2,bdry1);
grad = grad(:,:,1);
lhs(N_int+1,N_int+2) = grad;

[~,grad] = helm1d.green(k,int_chnkr,bdry2);
grad = grad(:,:,1);
lhs(N_int+2,1:N_int) = grad.*int_chnkr.wts(:).';
lhs(N_int+2,N_int+2) = 1/2;

[~,grad] = helm1d.green(k,bdry1,bdry2);
grad = grad(:,:,1);
lhs(N_int+2,N_int+1) = grad;

[~,grad2] = helm1d.green(k,src1,bdry1);
[~,grad3] = helm1d.green(k,src1,bdry2);
grad2 = grad2(:,:,1);
grad3 = grad3(:,:,1);
rhs = [skern(src1,int_chnkr) ; -grad2; -grad3];

mu = lhs \ rhs;
mu1 = mu(end-1);
mu2 = mu(end);
mu = mu(1:N_int);

G = chunkermat(tot_chnkr,gkern, opts);
usol = G(:,N_left+1:N_left+N_int)*mu + gkern(src1,tot_chnkr) + mu1*gkern(bdry1,tot_chnkr) + mu2*gkern(bdry2,tot_chnkr);


figure(10)
plot(int_chnkr.r(1,:),real(mu),int_chnkr.r(1,:),imag(mu))
title('\mu (interior density)')


figure(11)
plot(tot_chnkr.r(1,:),real(usol),tot_chnkr.r(1,:),imag(usol))
title('u (solution)')

% figure(12)
% uerrfk = usol - usol_fk;
% plot(tot_chnkr.r(1,:),real(uerrfk),tot_chnkr.r(1,:),imag(uerrfk))
% title('error')

% checking the solution

dmat = lege.dermat(16);
uext = usol(1:N_left);
u = squeeze(reshape(uext,size(left_chnkr.r(1,:,:))));
dd = squeeze(left_chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

f_recon = d2udx2(:)+k^2*u(:);

figure(12)
plot(left_chnkr.r(1,:),real(f_recon),left_chnkr.r(1,:),imag(f_recon))
title('Lu (= f)')

S = chunkermat(tot_chnkr,skern, opts2);
Su = S*usol;
int_err = usol(N_left+1:N_left+N_int) - Su(N_left+1:N_left+N_int);
figure(13);
plot(int_chnkr.r(1,:), real(int_err), int_chnkr.r(1,:), imag(int_err))
title('Interior error (u - S[u])')

figure(14);
umSu = squeeze(reshape(int_err,size(int_chnkr.r(1,:,:))));
dd = squeeze(int_chnkr.d(1,:,:));
dudx = dmat*umSu./dd;
d2udx2 = dmat*dudx./dd;
f_recon = d2udx2(:)+k^2*umSu(:);
plot(int_chnkr.r(1,:),real(f_recon),int_chnkr.r(1,:),imag(f_recon))
title('L(u-Su)')


% check that the BC is satisfied by the solution

u1 = u(:,end);
ds = left_chnkr.d(1,:,end).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(1,15);
dudn = pols.'*v2c*dmat*(u1./ds);


figure(15);
mu_err = mu - mu_true;
plot(int_chnkr.r(1,:),real(mu_err),int_chnkr.r(1,:),imag(mu_err))

%}

%% Method 4, use different representation 

lhs = zeros(N_int+2);

S_tot = chunkermat(tot_chnkr,skern,opts2);
S_int = chunkermat(int_chnkr,skern,opts2);
G_int = chunkermat(int_chnkr,gkern,opts);
G_tot = chunkermat(tot_chnkr,gkern,opts);

lhs(1:N_int,1:N_int) = eye(N_int) - k^2*G_int - S_int + k^2*S_tot(N_left+1:N_left+N_int,:)*G_tot(:,N_left+1:N_left+N_int);

lhs(1:N_int,N_int+1) = gkern(bdry1,int_chnkr) - S_tot(N_left+1:N_left+N_int,:)*gkern(bdry1,tot_chnkr);
lhs(1:N_int,N_int+2) = gkern(bdry2,int_chnkr) - S_tot(N_left+1:N_left+N_int,:)*gkern(bdry2,tot_chnkr);

[~,grad] = gkern(int_chnkr,bdry1);
grad = grad(:,:,1);
lhs(N_int+1,1:N_int) = -k^2*grad.*int_chnkr.wts(:).' ;
lhs(N_int+1,N_int+1) = -1/2;

[~,grad] = helm1d.green(k,bdry2,bdry1);
grad = grad(:,:,1);
lhs(N_int+1,N_int+2) = grad;

[~,grad] = gkern(int_chnkr,bdry2);
grad = grad(:,:,1);
lhs(N_int+2,1:N_int) = -k^2*grad.*int_chnkr.wts(:).' ;
lhs(N_int+2,N_int+2) = 1/2;

[~,grad] = gkern(bdry1,bdry2);
grad = grad(:,:,1);
lhs(N_int+2,N_int+1) = grad;

[~,grad2] = helm1d.green(k,src1,bdry1);
[~,grad3] = helm1d.green(k,src1,bdry2);

grad2 = grad2(:,:,1);
grad3 = grad3(:,:,1);

[~,grad2] = helm1d.green(k,left_chnkr,bdry1);
grad2 = grad2(:,:,1);

[~,grad3] = helm1d.green(k,left_chnkr,bdry2);
grad3 = grad3(:,:,1);

wts = left_chnkr.wts;

rhs = [- G(N_left+1:N_left+N_int,1:N_left)*f + S(N_left+1:N_left+N_int,:)*G(:,1:N_left)*f ; -grad2*(f.*wts(:)); -grad3*(f.*wts(:))]; 

mu = lhs \ rhs;
mu1 = mu(end-1);
mu2 = mu(end);
mu = mu(1:N_int);

figure(16)
plot(int_chnkr.r(1,:),real(mu),int_chnkr.r(1,:),imag(mu))
title('\mu (interior density)')

mu = [zeros(N_left,1); mu; zeros(N_right,1) ];
usol = mu - k^2*G_tot*mu + G(:,1:N_left)*f + mu1*gkern(bdry1,tot_chnkr) + mu2*gkern(bdry2,tot_chnkr);

figure(17)
plot(tot_chnkr.r(1,:),real(usol),tot_chnkr.r(1,:),imag(usol))
title('u (solution)')

% checking the solution

dmat = lege.dermat(16);
uext = usol(1:N_left);
u = squeeze(reshape(uext,size(left_chnkr.r(1,:,:))));
dd = squeeze(left_chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

f_recon = d2udx2(:)+k^2*u(:);

figure(18)
plot(left_chnkr.r(1,:),real(f-f_recon),left_chnkr.r(1,:),imag(f-f_recon))
title('Exterior error (Lu - f)')

Su = S_tot*usol;
int_err = usol(N_left+1:N_left+N_int) - Su(N_left+1:N_left+N_int);

figure(19);
plot(int_chnkr.r(1,:), real(int_err), int_chnkr.r(1,:), imag(int_err))
title('Interior error (u - S[u])')

% figure(20);
% umSu = squeeze(reshape(int_err,size(int_chnkr.r(1,:,:))));
% dd = squeeze(int_chnkr.d(1,:,:));
% dudx = dmat*umSu./dd;
% d2udx2 = dmat*dudx./dd;
% f_recon = d2udx2(:)+k^2*umSu(:);
% plot(int_chnkr.r(1,:),real(f_recon),int_chnkr.r(1,:),imag(f_recon))
% title('L(u-Su)')

% check that the BC is satisfied by the solution

u1 = u(:,end);
ds = left_chnkr.d(1,:,end).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(1,15);
dudn = pols.'*v2c*dmat*(u1./ds);

% figure(15);
% mu_err = mu - mu_true;
% plot(int_chnkr.r(1,:),real(mu_err),int_chnkr.r(1,:),imag(mu_err))

figure(21)
plot(tot_chnkr.r(1,:),real(usol - usol_analytic),tot_chnkr.r(1,:),imag(usol - usol_analytic))



