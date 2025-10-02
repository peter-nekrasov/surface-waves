%% Test S kernel using finite differences -- also ignore

L = 500;

n = 0;

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = L/10;
cparams.ta = 0;
cparams.tb = L; 

fcurve = @(t) [t(:).'; 0*t(:).']; 

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

figure(1)
plot(chnkr,'.')
title('chnkr')

sigma = 10;
f = exp(-(chnkr.r(1,:)-L/2).^2/(2*sigma^2)); 
f = f.';

figure(2)
plot(chnkr.r(1,:),real(f))
title('Original f')

skern =  @(s,t) axissymlap2d.kern(s,t,[0;0],'s'); 

opts2 = [];
opts2.sing = 'log';

S = chunkermat(chnkr,skern, opts2);

usol = -S*f;

figure(3)
plot(chnkr.r(1,:),real(usol),chnkr.r(1,:),imag(usol))
title('solution u')

dmat = lege.dermat(16);
u = squeeze(reshape(usol,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

rs = chnkr.r(1,:);
f1 = d2udx2(:) + 1./rs(:).*dudx(:) - n^2./rs(:).^2.*usol;

figure(21)
plot(chnkr.r(1,:),real(f1))
title('Reconstructed f')

tot_err = abs(f1 - f);

figure(22)
plot(chnkr.r(1,:),log10(tot_err))
title('Residual (interior and exterior)')
ylabel('log_{10} error')
hold on
