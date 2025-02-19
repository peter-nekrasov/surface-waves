% testing that the GF integrates to 0

L = 10000;

alpha = 2;
beta = 2+0.1i;
gamma = 3;
[rts,ejs] = helm1d.find_roots(alpha,beta,gamma);
k = rts(abs(angle(rts)) == min(abs(angle(rts))));

% create domains in chunkie 
cparams = [];
cparams.ifclosed = false;
cparams.maxchunklen = 5;
cparams.ta = 0;
cparams.tb = L; 

fcurve = @(t) [t(:).'; 0*t(:).']; 

chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

wts = chnkr.wts(:);

src = [];
src.r = [0;0];

[val] = helm1d.gshelm(rts,ejs,src,chnkr);

%[val] = helm1d.green(beta,src,chnkr);


int = val.'*wts