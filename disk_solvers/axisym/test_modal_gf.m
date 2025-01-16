%% Checking that the modal green's function satisfies the radial ODE

rs = lege.exps(16) + 3;
k = 2;

src = [];
src.r = [5; 0];

targ = [];
targ.r = [rs.'; 0*rs.'];

g = flex2d.modal_gf(k, src, targ, 0);
% g = flex2d.hkdiffgreen(k,src.r,targ.r);
dmat = lege.dermat(16);
dgdr = dmat*g;
d2gdr2 = dmat*dgdr;
d3gdr3 = dmat*d2gdr2;
d4gdr4 = dmat*d3gdr3;

err = d4gdr4+2./rs.*d3gdr3-1./rs.^2.*d2gdr2+1./rs.^3.*dgdr-k^4*g;

%% Checking that \int G f gives the correct solution to the ODE

k = 2+0.5i;

cparams = [];
cparams.ifclosed = false;
cparams.ta = 0;
cparams.tb = R;

fcurve = @(t) [t(:).';0*t(:).'];
cparams.maxchunklen = 1 / abs(k);
cparams.ta = 0;
cparams.tb = 10; % 10E5
chnkr = chunkerfunc(fcurve,cparams);
chnkr = sort(chnkr);

f = (chnkr.r(1,:)-5).'.*normpdf(5*chnkr.r(1,:)-25).';

figure(1);
t = tiledlayout(1,3);
nexttile 
plot(chnkr.r(1,:),f(:))
title('f(r)')

gkern =  @(s,t) flex2d.modal_gf(k, s, t, 0); 

opts = [];
opts.sing = 'log';
G = chunkermat(chnkr,gkern, opts);

sol = G*f;
nexttile 
plot(chnkr.r(1,:),real(sol),chnkr.r(1,:),imag(sol))
title('u = \int G f')
legend('real part','imaginary part')

dmat = lege.dermat(16);
sol = squeeze(reshape(sol,size(chnkr.r(1,:,:))));
dd = squeeze(chnkr.d(1,:,:));

u = sol;
dudr = dmat*u./dd;
d2udr2 = dmat*dudr./dd;
d3udr3 = dmat*d2udr2./dd;
d4udr4 = dmat*d3udr3./dd;

%hold on
%plot(chnkr.r(1,:),real(dudr(:)),chnkr.r(1,:),imag(dudr(:)))

rs = chnkr.r(1,:);
f_recon = d4udr4(:)+2./rs(:).*d3udr3(:)-1./rs(:).^2.*d2udr2(:)+1./rs(:).^3.*dudr(:)-k^4*u(:);

nexttile
plot(rs(50:end),real(f_recon(50:end)),rs(50:end),imag(f_recon(50:end)))
title('L u (= f)')




