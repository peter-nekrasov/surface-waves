L = 0;
k = 1+0.2i;

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 4 / abs(k);
cparams.ta = -100;
cparams.tb = -L; 

fcurve = @(t) [t(:).'; 0*t(:).'];

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

figure(1)
plot(chnkr,'.')
title('chnkr')


s = 2;
f = -(chnkr.r(1,:)+30).*exp(-(chnkr.r(1,:)+30).^2/(2*s^2)); % + vf(1)*besselh(n,k*chnkr.r(1,:)) + vf(2)*besselk(n,k*chnkr.r(1,:));
f = f.';

figure(2)
plot(chnkr.r(1,:), f)
title('f (RHS)')


targ = [];
targ.r = [-L; 0];

[val,grad] = helm1d.green(k,chnkr,targ);

wts = chnkr.wts;

grad = grad(:,:,1);
mu1 = 2*grad*(f.*wts(:));

% evaluation

gkern =  @(s,t) helm1d.green(k,s,t); 

opts = [];
opts.sing = 'removable';

G = chunkermat(chnkr,gkern, opts);

usol = G*f + mu1*val.';

figure(3)
plot(chnkr.r(1,:),real(usol),chnkr.r(1,:),imag(usol))
title('u (solution)')

% src = [];
% src.r = [-30;0];
% usol = helm1d.green(k,src,chnkr);

% checking the solution

dmat = lege.dermat(16);
u = squeeze(reshape(usol,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = chnkr.r(1,:);
f_recon = d2udx2(:)+k^2*u(:);

figure(4)
plot(rs,real(f_recon),rs,imag(f_recon))
title('Lu (= f)')

err = abs(f_recon(:) - f(:))/max(abs(f(:))).*sqrt(wts(:));

figure(5)
plot(rs,err)
title('relative error (residual)')

% check that the BC is satisfied by the solution

u1 = u(:,end);
ds = chnkr.d(1,:,end).';
[~,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(1,15);
dudn = pols.'*v2c*dmat*(u1./ds);