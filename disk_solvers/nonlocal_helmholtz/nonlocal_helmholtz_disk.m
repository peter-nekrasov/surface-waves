close all

L = 10; % width of polynya

n = 0; % mode number

alpha = 0.5;
beta = 1+0.2i;
gamma = 2;
[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
k = rts(abs(angle(rts)) == min(abs(angle(rts))));

right_bd = L;

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 8 / abs(k);
cparams.ta = right_bd;
cparams.tb = 80; 

fcurve = @(t) [t(:).'; 0*t(:).'];

right_chnkr = chunkerfunc(fcurve,cparams);
right_chnkr = sort(right_chnkr);

cparams.ta = 0;
cparams.tb = right_bd; 

cparams.maxchunklen = 8 / abs(k);
int_chnkr = chunkerfunc(fcurve,cparams);
int_chnkr = sort(int_chnkr);

tot_chnkr = merge([int_chnkr right_chnkr]);
tot_chnkr = sort(tot_chnkr);

N_int = int_chnkr.npt;
N_right = right_chnkr.npt;

figure(1)
plot(int_chnkr,'.')
hold on 
plot(right_chnkr,'.')
title('chnkr') 

L0 = 4*L; 

s = 1;
f = exp(-(right_chnkr.r(1,:)-L0).^2/(2*s^2));
f = f.';

figure(2)
plot(right_chnkr.r(1,:), f)
title('f (RHS)')

bdry2 = [];
bdry2.r = [right_bd; 0];

gskern =  @(s,t) helm2d.gsaxisym(rts,ejs,n,s,t); 
gphikern =  @(s,t) helm2d.gphiaxisym(rts,ejs,n,s,t); 

skern =  @(s,t) axissymlap2d.green(s.r,t.r,[0;0]); 

delgphikern = @(s,t) delphikern(rts,ejs,n,s,t);

opts = [];
opts.sing = 'removable';

opts2 = [];
opts2.sing = 'log';

lhs = zeros(N_int+1);

%S_tot = chunkermat(tot_chnkr,skern,opts2);
S_int = chunkermat(int_chnkr,skern,opts2);

Gs_int = chunkermat(int_chnkr,gskern,opts);
Gs_tot = chunkermat(tot_chnkr,gskern,opts);

return

Gphi_int = chunkermat(int_chnkr,gphikern,opts);
Gphi_tot = chunkermat(tot_chnkr,gphikern,opts);
delGphi_int = chunkermat(int_chnkr,delgphikern,opts2);
delGphi_tot = chunkermat(tot_chnkr,delgphikern,opts2);
SG_phi_int = -1/gamma*(S_int + beta/2*Gphi_int + alpha/2*delGphi_int);

lhs(1:N_int,1:N_int) = eye(N_int)/2 + beta/4*Gs_int + gamma/2*Gphi_int ...
    + gamma*S_int + gamma*beta/2*Gphi_int + gamma^2*SG_phi_int;

lhs(1:N_int,N_int+1) = 1/2*gskern(bdry1,int_chnkr) + gamma*gphikern(bdry1,int_chnkr);
lhs(1:N_int,N_int+2) = 1/2*gskern(bdry2,int_chnkr) + gamma*gphikern(bdry2,int_chnkr);

[~,gradgs] = gskern(int_chnkr,bdry1);
[~,gradgphi] = gphikern(int_chnkr,bdry1);
gradgs = gradgs(:,:,1);
gradgphi = gradgphi(:,:,1);
lhs(N_int+1,1:N_int) = beta/2*gradgs.*int_chnkr.wts(:).'+gamma*gradgphi.*int_chnkr.wts(:).' ;
lhs(N_int+1,N_int+1) = 1/alpha;
[~,grad] = gskern(bdry2,bdry1);
grad = grad(:,:,1);
lhs(N_int+1,N_int+2) = grad;

[~,gradgs] = gskern(int_chnkr,bdry2);
[~,gradgphi] = gphikern(int_chnkr,bdry2);
gradgs = gradgs(:,:,1);
gradgphi = gradgphi(:,:,1);
lhs(N_int+2,1:N_int) = beta/2*gradgs.*int_chnkr.wts(:).'+gamma*gradgphi.*int_chnkr.wts(:).' ;
lhs(N_int+2,N_int+2) = -1/alpha;

[~,grad] = gskern(bdry1,bdry2);
grad = grad(:,:,1);
lhs(N_int+2,N_int+1) = grad;

[~,grad2] = gskern(left_chnkr,bdry1);
grad2 = grad2(:,:,1);

[~,grad3] = gskern(left_chnkr,bdry2);
grad3 = grad3(:,:,1);

wts = left_chnkr.wts;

rhs = [1/2*Gs_tot(N_left+1:N_left+N_int,1:N_left)*f + gamma*Gphi_tot(N_left+1:N_left+N_int,1:N_left)*f; grad2*(f.*wts(:)); grad3*(f.*wts(:))]; 

mu = lhs \ rhs;
mu1 = mu(end-1);
mu2 = mu(end);
mu = mu(1:N_int);

%%

figure(16)
plot(int_chnkr.r(1,:),real(mu),int_chnkr.r(1,:),imag(mu))
title('\mu (interior density)')

mu = [zeros(N_left,1); mu; zeros(N_right,1) ];
usol = mu + beta/2*Gs_tot*mu + gamma*Gphi_tot*mu - Gs_tot(:,1:N_left)*f+ mu1*gskern(bdry1,tot_chnkr) + mu2*gskern(bdry2,tot_chnkr);

opts = [];
opts.sing = 'log';

Susol = -alpha/2*delGphi_tot*mu - Gphi_tot(:,1:N_left)*f+ mu1*gphikern(bdry1,tot_chnkr) + mu2*gphikern(bdry2,tot_chnkr);

figure(17)
plot(tot_chnkr.r(1,:),real(usol),tot_chnkr.r(1,:),imag(usol))
title('u (solution)')

% checking the solution

dmat = lege.dermat(16);
uleft = usol(1:N_left);
u = squeeze(reshape(uleft,size(left_chnkr.r(1,:,:))));
dd = squeeze(left_chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

left_err = 0.5*(alpha*d2udx2(:)+beta*u(:))+gamma*Susol(1:N_left) - f;
left_err = left_err ./ max(abs(uleft));

% figure(18)
% plot(left_chnkr.r(1,:),real(left_err),left_chnkr.r(1,:),imag(left_err))
% title('Left error (Lu - f)')

uint = usol(N_left+1:N_left+N_int);
int_err = 1/2*uint + gamma*Susol(N_left+1:N_left+N_int);
int_err = int_err ./ max(abs(uint));

% figure(20);
% plot(int_chnkr.r(1,:), real(int_err), int_chnkr.r(1,:), imag(int_err))
% title('Interior error 2 (u - S[u])')

uright = usol(N_left+N_int+1:end);
u = squeeze(reshape(uright,size(right_chnkr.r(1,:,:))));
dd = squeeze(right_chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

right_err = 0.5*(alpha*d2udx2(:)+beta*u(:))+gamma*Susol(N_left+N_int+1:end);
right_err = right_err ./ max(abs(uright));

% figure(21)
% plot(right_chnkr.r(1,:),real(right_err),right_chnkr.r(1,:),imag(right_err))
% title('Right error (Lu - f)')

%% 

tot_err = [abs(left_err); abs(int_err); abs(right_err)];
tot_err = tot_err.*sqrt(tot_chnkr.wts(:));


figure(22)

plot(tot_chnkr.r(1,:),log10(tot_err))
title('Residual (interior and exterior)')
ylabel('log_{10} error')
hold on

% figure(20);
% umSu = squeeze(reshape(int_err,size(int_chnkr.r(1,:,:))));
% dd = squeeze(int_chnkr.d(1,:,:));
% dudx = dmat*umSu./dd;
% d2udx2 = dmat*dudx./dd;
% f_recon = d2udx2(:)+k^2*umSu(:);
% plot(int_chnkr.r(1,:),real(f_recon),int_chnkr.r(1,:),imag(f_recon))
% title('L(u-Su)')

% check that the BC is satisfied by the solution

uleft = usol(1:N_left);
u = squeeze(reshape(uleft,size(left_chnkr.r(1,:,:))));
u1 = u(:,end);
ds = left_chnkr.d(1,:,end).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(1,15);
dudn1 = pols.'*v2c*dmat*(u1./ds);

% figure(15);
% mu_err = mu - mu_true;
% plot(int_chnkr.r(1,:),real(mu_err),int_chnkr.r(1,:),imag(mu_err))

uright = usol(end-N_right+1:end);
u = squeeze(reshape(uright,size(right_chnkr.r(1,:,:))));
u1 = u(:,1);
ds = left_chnkr.d(1,:,end).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);
dudn2 = pols.'*v2c*dmat*(u1./ds);


function hess = delphikern(rts,ejs,n,s,t)

    [~,~,hess] = helm2d.gphiaxisym(rts,ejs,n,s,t);
    
end