%% Checking (\Delta^2 + k^2) S[G(x,x_0)] = 1/(2*pi) log(|x-x_0|)

k = 1+0.2i;
x_0 = 5;

src = [];
src.r = [x_0;0];

cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 4 / abs(k);
cparams.ta = -100;
cparams.tb = 100; 

fcurve = @(t) [t(:).'; 0*t(:).'];

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);
N = chnkr.npt;

gkern =  @(s,t) helm1d.green(k,s,t); 
skern =  @(s,t) lap2d.kern(s,t,'s'); 

opts = [];
opts.sing = 'log';

G = chunkermat(chnkr,gkern, opts);
S = chunkermat(chnkr,skern, opts);
g = gkern(src,chnkr);
s = skern(src,chnkr);

t1 = S*g;

dmat = lege.dermat(16);
u = squeeze(reshape(t1,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));
dudx = dmat*u./dd;
d2udx2 = dmat*dudx./dd;

err = d2udx2(:) + k^2*u(:) - s(:);

plot(chnkr.r(1,:),d2udx2(:) + k^2*u(:),chnkr.r(1,:), s(:))
%plot(chnkr.r(1,1:N/2),real(err(1:N/2)), chnkr.r(1,1:N/2),imag(err(1:N/2)))

