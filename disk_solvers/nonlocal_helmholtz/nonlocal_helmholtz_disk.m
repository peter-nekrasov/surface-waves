close all

L = 10; % width of polynya

n = 1; % mode number (0 or 1)

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
cparams.tb = 60; 

fcurve = @(t) [t(:).'; 0*t(:).'];

ext_chnkr = chunkerfunc(fcurve,cparams);
ext_chnkr = sort(ext_chnkr);

cparams.ta = 0;
cparams.tb = right_bd; 

cparams.maxchunklen = 8 / abs(k);
int_chnkr = chunkerfunc(fcurve,cparams);
int_chnkr = sort(int_chnkr);

tot_chnkr = merge([int_chnkr ext_chnkr]);
tot_chnkr = sort(tot_chnkr);

N_int = int_chnkr.npt;
N_ext = ext_chnkr.npt;

figure(1)
plot(int_chnkr,'.')
hold on 
plot(ext_chnkr,'.')
title('chnkr') 

bdry = [];
bdry.r = [right_bd; 0];

gskern =  @(s,t) helm2d.gsaxisym(rts,ejs,n,s,t); 
gphikern =  @(s,t) helm2d.gphiaxisym(rts,ejs,n,s,t); 
skern =  @(s,t) axissymlap2d.green(n,s.r,t.r,[0;0]); 
delgphikern = @(s,t) delphikern(rts,ejs,n,s,t);
gradkern =  @(s,t) gradgskern(rts,ejs,n,s,t); 

opts = [];
opts.sing = 'removable';

opts2 = [];
opts2.sing = 'log';

eval_opts = [];
eval_opts.forcesmooth = true;

%S_tot = chunkermat(tot_chnkr,skern,opts)/(4*pi^2);
S_i2i = chunkermat(int_chnkr,skern,opts)/(4*pi^2);
Gs_i2i = chunkermat(int_chnkr,gskern,opts);
Gs_t2t = chunkermat(tot_chnkr,gskern,opts);
Gphi_i2i = chunkermat(int_chnkr,gphikern,opts);
Gphi_t2t = chunkermat(tot_chnkr,gphikern,opts);
delGphi_i2i = chunkermat(int_chnkr,delgphikern,opts2);
delGphi_t2t = chunkermat(tot_chnkr,delgphikern,opts2);

Gs_e2i = chunkerkernevalmat(ext_chnkr,gskern,int_chnkr,eval_opts);
Gphi_e2i = chunkerkernevalmat(ext_chnkr,gphikern,int_chnkr,eval_opts);

%%

L0 = 2*L; 

s = 1;
f = exp(-(ext_chnkr.r(1,:)-L0).^2/(2*s^2));
f = f.';

figure(2)
plot(ext_chnkr.r(1,:), f)
title('f (RHS)')

lhs = zeros(N_int+1);

SG_phi_i2i = -1/gamma*(S_i2i + beta/2*Gphi_i2i + alpha/2*delGphi_i2i);

lhs(1:N_int,1:N_int) = eye(N_int)/2 + beta/4*Gs_i2i + gamma/2*Gphi_i2i ...
    + gamma*S_i2i + gamma*beta/2*Gphi_i2i + gamma^2*SG_phi_i2i;


lhs(1:N_int,N_int+1) = 1/2*gskern(bdry,int_chnkr) + gamma*gphikern(bdry,int_chnkr);

%lhs(1:N_int,N_int+2) = 1/2*gskern(bdry,int_chnkr) + gamma*gphikern(bdry,int_chnkr);
%
% [~,gradgs] = gskern(int_chnkr,bdry1);
% [~,gradgphi] = gphikern(int_chnkr,bdry1);
% gradgs = gradgs(:,:,1);
% gradgphi = gradgphi(:,:,1);
% lhs(N_int+1,1:N_int) = beta/2*gradgs.*int_chnkr.wts(:).'+gamma*gradgphi.*int_chnkr.wts(:).' ;
% lhs(N_int+1,N_int+1) = 1/alpha;
% [~,grad] = gskern(bdry,bdry1);
% grad = grad(:,:,1);
% lhs(N_int+1,N_int+2) = grad;
% 

dmat = lege.dermat(16);
ds = ext_chnkr.d(1,:,1).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);

% [~,gradgs] = gskern(int_chnkr,bdry);
% [~,gradgphi] = gphikern(int_chnkr,bdry);
lhs(N_int+1,1:N_int) = (pols.'*v2c*dmat./ds.')*(beta/2*Gs_t2t(N_int+1:N_int+16,1:N_int)+gamma*Gphi_t2t(N_int+1:N_int+16,1:N_int)); % beta/2*gradgs.*int_chnkr.wts(:).'+gamma*gradgphi.*int_chnkr.wts(:).' ;

bterm = gskern(bdry,ext_chnkr);
bterm = bterm(1:16);
lhs(N_int+1,N_int+1) = pols.'*v2c*dmat*(bterm(:)./ds);

grad = gradkern(bdry,bdry);

lhs(N_int+1,N_int+1) = -1/alpha + grad;

% [~,grad] = gskern(bdry1,bdry);
% grad = grad(:,:,1);
% lhs(N_int+2,N_int+1) = grad;
% 
% [~,grad2] = gskern(left_chnkr,bdry1);
% grad2 = grad2(:,:,1);
% 

%grad2 = chunkerkerneval(ext_chnkr,gradkern,bdry,eval_opts);

wts = ext_chnkr.wts;
rhs = [1/2*Gs_e2i*f + gamma*Gphi_e2i*f; (pols.'*v2c*dmat./ds.')*Gs_t2t(N_int+1:N_int+16,N_int+1:end)*f]; %grad2*(f.*wts(:))]; %; grad3*(f.*wts(:))]; 

mu = lhs \ rhs;
mu1 = mu(end);
mu = mu(1:N_int);

figure(16)
plot(int_chnkr.r(1,:),real(mu),int_chnkr.r(1,:),imag(mu))
title('\mu (interior density)')

mu = [mu; zeros(N_ext,1)];
usol = mu + beta/2*Gs_t2t*mu + gamma*Gphi_t2t*mu - Gs_t2t(:,N_int+1:end)*f + mu1*gskern(bdry,tot_chnkr);

Susol = -alpha/2*delGphi_t2t*mu - Gphi_t2t(:,N_int+1:end)*f + mu1*gphikern(bdry,tot_chnkr); % + mu2*gphikern(bdry,tot_chnkr);

figure(17)
plot(tot_chnkr.r(1,:),real(usol),tot_chnkr.r(1,:),imag(usol))
title('u (solution)')


% checking the solution


uint = usol(1:N_int);
int_err = 1/2*uint + gamma*Susol(1:N_int);
int_err = int_err ./ max(abs(uint));

figure(20);
plot(int_chnkr.r(1,:), real(int_err), int_chnkr.r(1,:), imag(int_err))
title('Interior error')

uright = usol(N_int+1:end);
u = squeeze(reshape(uright,size(ext_chnkr.r(1,:,:))));
dd = squeeze(ext_chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = ext_chnkr.r(1,:);
f1 = alpha/2*d2udx2(:) + alpha/2./rs(:).*dudx(:) - alpha/2*n^2./rs(:).^2.*u(:) + beta/2*u(:) + gamma*Susol(N_int+1:end);

ext_err = abs(f-f1);

% figure(21)
% plot(right_chnkr.r(1,:),real(right_err),right_chnkr.r(1,:),imag(right_err))
% title('Right error (Lu - f)')


tot_err = [abs(int_err); abs(ext_err)];
tot_err = tot_err.*sqrt(tot_chnkr.wts(:));


figure(21)

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

% figure(15);
% mu_err = mu - mu_true;
% plot(int_chnkr.r(1,:),real(mu_err),int_chnkr.r(1,:),imag(mu_err))

uright = usol(N_int+1:end);
u = squeeze(reshape(uright,size(ext_chnkr.r(1,:,:))));
u1 = u(:,1);
ds = ext_chnkr.d(1,:,1).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);
dudn = pols.'*v2c*dmat*(u1./ds)

%%

cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 8 / abs(k);
cparams.ta = 1.2*right_bd;
cparams.tb = 60; 

ext_chnkr = chunkerfunc(fcurve,cparams);
ext_chnkr = sort(ext_chnkr);

cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 2 / abs(k);
cparams.ta = right_bd;
cparams.tb = 1.2*right_bd;

ref_chnkr = chunkerfunc(fcurve,cparams);
ref_chnkr = sort(ref_chnkr);

Gs_e2r = chunkerkernevalmat(ext_chnkr,gskern,ref_chnkr,eval_opts);
Gphi_e2r = chunkerkernevalmat(ext_chnkr,gphikern,ref_chnkr,eval_opts);

Gs_i2r = chunkerkernevalmat(int_chnkr,gskern,ref_chnkr,eval_opts);
Gphi_i2r = chunkerkernevalmat(int_chnkr,gphikern,ref_chnkr,eval_opts);
delGphi_i2r = chunkerkernevalmat(int_chnkr,delgphikern,ref_chnkr,eval_opts);

%%

mu_int = mu(1:N_int);
usol = beta/2*Gs_i2r*mu_int + gamma*Gphi_i2r*mu_int - Gs_e2r*f + mu1*gskern(bdry,ref_chnkr);
Susol = -alpha/2*delGphi_i2r*mu_int - Gphi_e2r*f + mu1*gphikern(bdry,ref_chnkr); 

dmat = lege.dermat(16);
u = squeeze(reshape(usol,size(ref_chnkr.r(1,:,:))));
dd = squeeze(ref_chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = ref_chnkr.r(1,:);
f1 = alpha/2*d2udx2(:) + alpha/2./rs(:).*dudx(:) - alpha/2*n^2./rs(:).^2.*u(:) + beta/2*u(:) + gamma*Susol(:);

f2 = exp(-(ref_chnkr.r(1,:)-L0).^2/(2*s^2));
ext_err = abs(f1-f2.').*sqrt(ref_chnkr.wts(:));

figure(23)
plot(ref_chnkr.r(1,:),real(usol(:)),ref_chnkr.r(1,:),imag(usol(:)))
title('Solution (u)')
hold on

figure(24)
plot(ref_chnkr.r(1,:),log10(ext_err),'x-')
title('Residual (closest panel)')
ylabel('log_{10} error')
hold on

u = squeeze(reshape(usol,size(ref_chnkr.r(1,:,:))));
u1 = u(:,1);
ds = ref_chnkr.d(1,:,1).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);
dudn = pols.'*v2c*dmat*(u1./ds);

load('gong.mat')
sound(y)


function hess = delphikern(rts,ejs,n,s,t)

    [~,~,hess] = helm2d.gphiaxisym(rts,ejs,n,s,t);
    
end

function grad = gradgskern(rts,ejs,n,s,t)

    [~,grad] = helm2d.gsaxisym(rts,ejs,n,s,t);
    
end