[x,w] = lege.exps(16);
f = x.^2;

dmat = lege.dermat(16);
dfdx = dmat*f;
d2fdx2 = dmat*dmat*f;

plot(x,f,x,dfdx,x,d2fdx2)

%% 

cparams = [];
cparams.ifclosed = false;
cparams.ta = 0;
cparams.tb = 1;
pref = [];
pref.nchmax = 1;

fcurve = @(t) [t(:).';0*t(:).'];
int_chnkr = chunkerfunc(fcurve,cparams,pref);

f = int_chnkr.r(1,:).^2;
dd = int_chnkr.d(1,:);

dmat = lege.dermat(16);
dfdx = dmat*f(:)./dd(:);
d2fdx2 = dmat*dfdx ./ dd(:);

plot(int_chnkr.r(1,:),0*int_chnkr.r(1,:),'x-')
hold on
plot(int_chnkr.r(1,:),f,int_chnkr.r(1,:),dfdx,int_chnkr.r(1,:),d2fdx2)

%% 

cparams = [];
cparams.maxchunklen = 8 / abs(k);
cparams.ifclosed = false;
cparams.ta = 0;
cparams.tb = sqrt(2);

fcurve = @(t) [t(:).';0*t(:).'];
int_chnkr = chunkerfunc(fcurve,cparams);

cparams.ta = sqrt(2);
cparams.tb = 5; % 10E5
ext_chnkr = chunkerfunc(fcurve,cparams);

chnkr = merge([int_chnkr ext_chnkr]);

dmat = repmat(lege.dermat(16),chnkr.nch,chnkr.nch);
f = chnkr.r(1,:).^2;
dd = chnkr.d(1,:);
dfdx = dmat*f(:)./d(:);
d2fdx2 = dmat*dmat*f(:)./d(:);

plot(chnkr.r(1,:),f,chnkr.r(1,:),dfdx,chnkr.r(1,:),d2fdx2)
