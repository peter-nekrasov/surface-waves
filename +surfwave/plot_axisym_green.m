%% Plotting axisymmetric kernels

close all

alpha = 0.5;
beta = 1.5+0.5i;
gamma = 2.5;
[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
k = rts(abs(angle(rts)) == min(abs(angle(rts))));

n = 0;

src = [];
src.r = [8;0];

h = 0.01;
targ = [];
targ.r = [0:h:20; 0*(0:h:20)];

tic
[valgs,gradgs,lapgs] = helm2d.gsaxisym(rts,ejs,n,src,targ);
toc

figure(1)
tiledlayout(1,2);
nexttile
plot(targ.r(1,:),real(valgs),targ.r(1,:),imag(valgs))
title('G_S (nonlocal Helmholtz)')

tic
[valgphi,gradgphi,lapgphi] = helm2d.gphiaxisym(rts,ejs,n,src,targ);
toc

nexttile
plot(targ.r(1,:),real(valgphi),targ.r(1,:),imag(valgphi))
title('G_\phi (or S[G_S])')

% 
% [val] = chnk.axissymhelm2d.green(k,src.r,targ.r,[0;0]);
% 
% nexttile
% plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
% title('G (Helmholtz)')

%%

figure(2)
tiledlayout(1,2);
nexttile
plot(targ.r(1,:),real(gradgs),targ.r(1,:),imag(gradgs))
title('\nabla G_S (nonlocal Helmholtz)')

nexttile
plot(targ.r(1,:),real(gradgphi),targ.r(1,:),imag(gradgphi))
title('\nabla G_\phi')

%%

figure(3)
tiledlayout(1,2);
nexttile
plot(targ.r(1,:),real(lapgs),targ.r(1,:),imag(lapgs))
title('\Delta G_S (nonlocal Helmholtz)')


nexttile
plot(targ.r(1,:),real(lapgphi),targ.r(1,:),imag(lapgphi))
title('\Delta G_\phi')
