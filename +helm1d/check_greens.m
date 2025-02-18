%% test 1 - finite difference

h = 0.01;

src = [];
src.r = [5;0];

targ = [];
targ.r = [0:h:8*h; (0:h:8*h)*0];

d2 = [-1/560 8/315 -1/5 8/5 -205/72 8/5 -1/5 8/315 -1/560]/h^2;

% test on helmholtz green function first

k = 2;
[g,~] = helm1d.green(k,src,targ);

err = d2*g + k^2*g(5); % passed


% now check nonlocal greens function

beta = 1;
gamma = 2;
[rts,ejs] = helm1d.find_roots(beta,gamma);

[gs,~] = helm1d.gshelm(rts,ejs,src,targ);
[gphi,~] = helm1d.gphihelm(rts,ejs,src,targ);

err = d2*gs + beta*gs(5) + gamma*gphi(5); % passed

%% test 2 - part 1 - solve (\Delta + k^2)u = f 
% and check consistency using spectral differentiation

L = 200;

beta = 1;
gamma = 2;
[rts,ejs] = helm1d.find_roots(beta,gamma);
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

gkern =  @(s,t) helm1d.green(k,s,t); 

opts = [];
opts.sing = 'removable';

G = chunkermat(chnkr,gkern, opts);
u = G*f;

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




%% test 2 - part 2 - solve (\Delta + \beta)u + \gamma S[u] = f 
% and check consistency using spectral differentiation

L = 800;

beta = 1+0.5i;
gamma = 1;
[rts,ejs] = helm1d.find_roots(beta,gamma);
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

sigma = 5;
f = exp(-(chnkr.r(1,:)-L/2).^2/(2*sigma^2)); 
f = f.';

skern =  @(s,t) chnk.lap2d.kern(s,t,'s'); 
gskern =  @(s,t) helm1d.gshelm(rts,ejs,s,t); 
gphikern =  @(s,t) helm1d.gphihelm(rts,ejs,s,t); 

opts = [];
opts.sing = 'removable';

opts2 = [];
opts2.sing = 'log';

Gs = chunkermat(chnkr,gskern, opts);
Gphi = chunkermat(chnkr,gphikern, opts);
S = chunkermat(chnkr,skern, opts2);

u = Gs*f;
Su = Gphi*f;
Su_prod = S*Gs*f;

figure(2)
plot(chnkr.r(1,:),real(Su),chnkr.r(1,:),imag(Su))
title('Su')

figure(3)
plot(chnkr.r(1,:),real(Su_prod),chnkr.r(1,:),imag(Su_prod))
title('S*u')

figure(4)
plot(chnkr.r(1,:),real(u),chnkr.r(1,:),imag(u))
title('u')

figure(5)
tiledlayout(1,3)
nexttile
plot(chnkr.r(1,:), f)
title('f (RHS)')

return 

dmat = lege.dermat(16);
u = squeeze(reshape(u,size(chnkr.r(1,:,:))));
Su = squeeze(reshape(Su,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = chnkr.r(1,:);
f_recon = d2udx2(:)+beta*u(:)+gamma*Su(:); %

nexttile
plot(rs,real(f_recon))
title('(\Delta + \beta)u + \gamma S[u] (= f_0)')

wts = chnkr.wts(:);
err = abs(f_recon(:) - f(:))/max(abs(f(:)));

nexttile
plot(rs,real(err))
title('relative error |f_0 - f|/max(|f|)')


%% test 3 - part 1 - (\Delta + k^2)u = f, but evaluate solution at discrete 
% set of points so that we can use finite difference instead of spectral diff
% INCOMPLETE

L = 200;

beta = 1;
gamma = 2;
[rts,ejs] = helm1d.find_roots(beta,gamma);
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

s = 2;
f = exp(-(chnkr.r(1,:)-L/3).^2/(2*s^2)); 
f = f.';

tiledlayout(1,3)
nexttile
plot(chnkr.r(1,:), f)
title('f (RHS)')

gkern =  @(s,t) helm1d.green(k,s,t); 


h = 0.01;

targ = [];
loc = 65;
targ.r = [loc:h:(loc+8*h); (loc:h:(loc+8*h))*0];

% wts = chnkr.wts(:);
% 
% Gs = gskern(chnkr,targ).*wts.';
% Gphi = gphikern(chnkr,targ).*wts.';

u = chunkerkerneval(chnkr,gkern,f,targ);

d2 = [-1/560 8/315 -1/5 8/5 -205/72 8/5 -1/5 8/315 -1/560]/h^2;

err = d2*u+k^2*u(5)



%% same thing, but for (\Delta + \beta)u + \gamma S[u] = f,

L = 200;

beta = 1;
gamma = 2;
[rts,ejs] = helm1d.find_roots(beta,gamma);
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

s = 2;
f = exp(-(chnkr.r(1,:)-L/3).^2/(2*s^2)); 
f = f.';

tiledlayout(1,3)
nexttile
plot(chnkr.r(1,:), f)
title('f (RHS)')

gskern =  @(s,t) helm1d.gshelm(rts,ejs,s,t); 
gphikern =  @(s,t) helm1d.gphihelm(rts,ejs,s,t); 

h = 0.01;

targ = [];
loc = 66;
targ.r = [loc:h:(loc+8*h); (loc:h:(loc+8*h))*0];

% wts = chnkr.wts(:);
% 
% Gs = gskern(chnkr,targ).*wts.';
% Gphi = gphikern(chnkr,targ).*wts.';

u = chunkerkerneval(chnkr,gskern,f,targ);
Su = chunkerkerneval(chnkr,gphikern,f,targ);

d2 = [-1/560 8/315 -1/5 8/5 -205/72 8/5 -1/5 8/315 -1/560]/h^2;

err = d2*u+beta*u(5)+gamma*Su(5)