%% Integrating the axisymmetric G against some rhs f
close all

n = 0; % mode number

alpha = 0.5;
beta = 1+0.2i;
gamma = 2;
[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
k = rts(abs(angle(rts)) == min(abs(angle(rts))));


% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 8 / abs(k);
cparams.ta = 0;
cparams.tb = 70; 

fcurve = @(t) [t(:).'; 0*t(:).'];

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

figure(1)
plot(chnkr,'.')
hold on 

s = 1;
f = exp(-(chnkr.r(1,:)-30).^2/(2*s^2));
f = f.';

figure(2)
plot(chnkr.r(1,:), f)
title('f (RHS)')

gskern =  @(s,t) helm2d.gsaxisym(rts,ejs,n,s,t); 
gphikern =  @(s,t) helm2d.gphiaxisym(rts,ejs,n,s,t); 

opts = [];
opts.sing = 'removable';

Gs = chunkermat(chnkr,gskern,opts);
Gphi = chunkermat(chnkr,gphikern,opts);

%%

usol = - Gs*f;
Susol = - Gphi*f;

figure(17)
plot(chnkr.r(1,:),real(usol),chnkr.r(1,:),imag(usol))
title('u (solution)')

% checking the solution

dmat = lege.dermat(16);
u = squeeze(reshape(usol,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = chnkr.r(1,:);
f1 = alpha/2*d2udx2(:) + alpha/2./rs(:).*dudx(:) - alpha/2*n^2./rs(:).^2.*usol + beta/2*usol + gamma*Susol;

figure(21)
plot(chnkr.r(1,:),real(f1))
title('Reconstructed f')

tot_err = abs(f1 - f);

figure(22)
plot(chnkr.r(1,:),log10(tot_err))
title('Residual (interior and exterior)')
ylabel('log_{10} error')
hold on