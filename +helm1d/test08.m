%% test 1 - finite difference

h = 0.01;

src = [];
src.r = [5;0];

targ = [];
targ.r = [0:h:8*h; (0:h:8*h)*0];

d2 = [-1/560 8/315 -1/5 8/5 -205/72 8/5 -1/5 8/315 -1/560]/h^2;

% test on helmholtz green function first

k = 2;
[g,~] = chnk.helm1d.green(k,src.r,targ.r);

err = d2*g + k^2*g(5); % passed


% now check nonlocal greens function

alpha = 1;
beta = 1;
gamma = 2;
[rts,ejs] = helm1d.find_roots(alpha,beta,gamma);

[gs,~] = helm1d.gshelm(rts,ejs,src,targ);
[gphi,~] = helm1d.gphihelm(rts,ejs,src,targ);

err = 1/2*d2*gs + 1/2*beta*gs(5) + gamma*gphi(5) % passed

%% test 2 - part 1 - solve (\Delta + k^2)u = f 
% and check consistency using spectral differentiation

L = 200;

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 8 / abs(k);
cparams.ta = 0;
cparams.tb = L; 

fcurve = @(t) [t(:).'; 0*t(:).']; 

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

figure(1)
plot(chnkr,'.')
title('chnkr')

s = 4;
f = exp(-(chnkr.r(1,:)-L/3).^2/(2*s^2)); 
f = f.';

gkern =  @(s,t) 1i*chnk.helm1d.green(k,s.r,t.r)/(k^2); 

opts = [];
opts.sing = 'removable';

G = chunkermat(chnkr,gkern, opts);
u = -G*f;

figure(1)
plot(chnkr.r(1,:),real(u),chnkr.r(1,:),imag(u))
title('u')

figure(2)
tiledlayout(1,3)
nexttile
plot(chnkr.r(1,:), f)
title('f (RHS)')

dmat = lege.dermat(16);
u = squeeze(reshape(u,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = chnkr.r(1,:);
f_recon = d2udx2(:)+k^2*u(:);

nexttile
plot(rs,real(f_recon))
title('(\Delta + \beta)u (= f_0)')

wts = chnkr.wts(:);
err = abs(f_recon(:) - f(:))/max(abs(f(:))).*sqrt(wts);

nexttile
plot(rs,real(err))
title('relative error |f_0 - f|/max(|f|)')




%% test 2 - part 2 - solve (\alpha \Delta + \beta)u + \gamma S[u] = f 
% and check consistency using spectral differentiation

L = 200;

alpha = 0.5;
beta = 1;
gamma = 2;
[rts,ejs] = helm1d.find_roots(alpha,beta,gamma);
k = rts(abs(angle(rts)) == min(abs(angle(rts))));

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 8 / abs(k);
cparams.ta = 0;
cparams.tb = L; 

fcurve = @(t) [t(:).'; 0*t(:).']; 

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

figure(1)
plot(chnkr,'.')
title('chnkr')

s = 4;
f = exp(-(chnkr.r(1,:)-L/3).^2/(2*s^2));  
f = f.';

gskern =  @(s,t) helm1d.gshelm(rts,ejs,s,t); 
gphikern =  @(s,t) helm1d.gphihelm(rts,ejs,s,t); 

opts = [];
opts.sing = 'removable';

Gs = chunkermat(chnkr,gskern, opts);
Gphi = chunkermat(chnkr,gphikern, opts);

opts = [];
opts.sing = 'removable';


usol = -Gs*f;
Su = -Gphi*f;

figure(2)
plot(chnkr.r(1,:),real(Su),chnkr.r(1,:),imag(Su))
title('Su')


figure(3)
plot(chnkr.r(1,:),real(usol),chnkr.r(1,:),imag(usol))
title('u')



dmat = lege.dermat(16);
u = squeeze(reshape(usol,size(chnkr.r(1,:,:))));
Su = squeeze(reshape(Su,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

figure(5)
plot(chnkr.r(1,:),real(d2udx2(:)),chnkr.r(1,:),imag(d2udx2(:)))
title('\Delta u')

rs = chnkr.r(1,:);
f_recon = 1/2*alpha*d2udx2(:)+1/2*beta*u(:)+gamma*Su(:); %

figure(6)
tiledlayout(1,3)
nexttile
plot(chnkr.r(1,:), f)
title('f (RHS)')

nexttile
plot(rs,real(f_recon))
title('1/2(\alpha\Delta + \beta)u + \gamma S[u] (= f_0)')

wts = chnkr.wts(:);
err = abs(f_recon(:) - f(:))/max(abs(f(:))).*sqrt(chnkr.wts(:));

nexttile
plot(rs,real(err))
title('relative error |f_0 - f|/max(|f|)')