%% Test S[G_S] = G_phi
% Expect 1/R convergence, where R is the length of the domain

L = 80;

n = 0;

alpha = 2;
beta = 3+1i;
gamma = 4;
[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
k = rts(abs(angle(rts)) == min(abs(angle(rts))));

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 16 / abs(k);
cparams.ta = 0;
cparams.tb = L; 

fcurve = @(t) [t(:).'; 0*t(:).']; 

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

figure(1)
plot(chnkr,'.')
title('chnkr')

sigma = 2;
f = exp(-(chnkr.r(1,:)-L/2).^2/(2*sigma^2)); 
f = f.';

skern =  @(s,t) axissymlap2d.green(s.r,t.r,'s'); 
gskern =  @(s,t) helm2d.gsaxisym(rts,ejs,n,s,t); 
gphikern =  @(s,t) helm2d.gphiaxisym(rts,ejs,n,s,t); 

opts = [];
opts.sing = 'removable';

opts2 = [];
opts2.sing = 'log';

Gs = chunkermat(chnkr,gskern, opts);
Gphi = chunkermat(chnkr,gphikern, opts);
S = chunkermat(chnkr,skern, opts2);

SG1 = -S*Gs*f/(2*pi)^2;
SG2 = -Gphi*f;

figure(2)
plot(chnkr.r(1,:),real(SG1),chnkr.r(1,:),imag(SG1))
title('S[G]')

figure(3)
plot(chnkr.r(1,:),real(SG2),chnkr.r(1,:),imag(SG2))
title('G_\phi')

figure(4)
err = abs(SG1 - SG2);
plot(chnkr.r(1,:),err)
title('Absolute error')

figure(5)
tiledlayout(1,3)
nexttile
plot(chnkr.r(1,:), f)
title('f (RHS)')




%% Test S[G_phi] kernel 
% See notes for identity

L = 80;

n = 0;

alpha = 2;
beta = 3+1i;
gamma = 4;
[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
k = rts(abs(angle(rts)) == min(abs(angle(rts))));

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 16 / abs(k);
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

skern =  @(s,t) axissymlap2d.green(s.r,t.r,'s'); 
gskern =  @(s,t) helm2d.gsaxisym(rts,ejs,n,s,t); 
gphikern =  @(s,t) helm2d.gphiaxisym(rts,ejs,n,s,t); 
delgphikern =  @(s,t) delgkern(rts,ejs,n,s,t); 

opts = [];
opts.sing = 'removable';

opts2 = [];
opts2.sing = 'log';

opts3 = [];
opts3.sing = 'log';

Gs = chunkermat(chnkr,gskern, opts);
Gphi = chunkermat(chnkr,gphikern, opts);
delGphi = chunkermat(chnkr,delgphikern, opts3);
S = chunkermat(chnkr,skern,opts2);

%%

SGphi1 = S*Gphi*f/(2*pi)^2;
SGphi2 = -1/gamma*(S/(2*pi)^2 + 0.5*beta*Gphi + 0.5*alpha*delGphi)*f;

figure(2)
plot(chnkr.r(1,:),real(SGphi1),chnkr.r(1,:),imag(SGphi1))
title('S[G_phi]')

figure(3)
plot(chnkr.r(1,:),real(SGphi2),chnkr.r(1,:),imag(SGphi2))
title('1/\gamma(S - \beta*G_\phi - \Delta G_\phi)')

figure(4)
err = abs(SGphi1 - SGphi2);
plot(chnkr.r(1,:),err)
title('Absolute error')

figure(5)
tiledlayout(1,3)
nexttile
plot(chnkr.r(1,:), f)
title('f (RHS)')

function hess = delgkern(rts,ejs,n,s,t)

    [~,~,hess] = helm2d.gphiaxisym(rts,ejs,n,s,t);
    
end