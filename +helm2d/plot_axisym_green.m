%% Plotting axisymmetric kernels

close all

alpha = 0.5;
beta = 1.5+0.5i;
gamma = 2.5;
[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
k = rts(abs(angle(rts)) == min(abs(angle(rts))));

n = 2;

src = [];
src.r = [8;0];

targ = [];
targ.r = [0:0.1:20; (0:0.1:20)*0];

[valgs,gradgs] = helm2d.gsaxisym(rts,ejs,n,src,targ);

figure(1)
tiledlayout(1,3);
nexttile
plot(targ.r(1,:),real(valgs),targ.r(1,:),imag(valgs))
title('G_S (nonlocal Helmholtz)')


[valgphi] = helm2d.gphiaxisym(rts,ejs,n,src,targ);

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
tiledlayout(1,3);
nexttile
plot(targ.r(1,:),real(gradgs),targ.r(1,:),imag(gradgs))
title('\nabla G_S (nonlocal Helmholtz)')

return

% [~,grad] = helm2d.gphihelm(rts,ejs,src,targ);
% 
% nexttile
% plot(targ.r(1,:),real(grad),targ.r(1,:),imag(grad))
% title('\nabla G_\phi (or S[G_S])')

[~,grad] = chnk.helm2d.green(k,src.r,targ.r);
grad = grad(:,:,1);

nexttile
plot(targ.r(1,:),real(grad),targ.r(1,:),imag(grad))
title('\nabla G (Helmholtz)')

%%

targ.r = [-1:0.001:1; (-1:0.001:1)*0];

[~,~,hess] = helm2d.gshelm(rts,ejs,src,targ);
hess = hess(:,:,1);

figure(3)
tiledlayout(1,3);
nexttile
plot(targ.r(1,:),real(hess),targ.r(1,:),imag(hess))
title('\Delta G_S (nonlocal Helmholtz)')

% [~,~,hess] = helm2d.gphihelm(rts,ejs,src,targ);
% 
% nexttile
% plot(targ.r(1,:),real(hess),targ.r(1,:),imag(hess))
% title('\Delta G_\phi (or S[G_S])')

[~,~,hess] = chnk.helm2d.green(k,src.r,targ.r);
hess = hess(:,:,1);

nexttile
plot(targ.r(1,:),real(hess),targ.r(1,:),imag(hess))
title('\Delta G (Helmholtz)')
