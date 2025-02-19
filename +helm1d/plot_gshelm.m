close all

alpha = 0.5;
beta = 1+0.5i;
gamma = 2;
[rts,ejs] = helm1d.find_roots(alpha,beta,gamma);

src = [];
src.r = [0;0];

targ = [];
targ.r = [-20:0.01:20; (-20:0.01:20)*0];

[val,~] = helm1d.gshelm(rts,ejs,src,targ);

figure(1)
tiledlayout(1,3);
nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G_S (nonlocal Helmholtz)')


[val,~] = helm1d.gphihelm(rts,ejs,src,targ);

nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G_\phi (or S[G_S])')

k = rts(abs(angle(rts)) == min(abs(angle(rts))));
[val,~] = chnk.helm1d.green(k,src.r,targ.r);
val = 1i*val(:,:,1)./(2*k);

nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G (Helmholtz)')

%%

[~,grad] = helm1d.gshelm(rts,ejs,src,targ);

figure(2)
tiledlayout(1,3);
nexttile
plot(targ.r(1,:),real(grad),targ.r(1,:),imag(grad))
title('\nabla G_S (nonlocal Helmholtz)')


[~,grad] = helm1d.gphihelm(rts,ejs,src,targ);

nexttile
plot(targ.r(1,:),real(grad),targ.r(1,:),imag(grad))
title('\nabla G_\phi (or S[G_S])')

k = rts(abs(angle(rts)) == min(abs(angle(rts))));
[~,grad] = chnk.helm1d.green(k,src.r,targ.r);
grad = 1i*grad(:,:,1)./(2*k);

nexttile
plot(targ.r(1,:),real(grad),targ.r(1,:),imag(grad))
title('\nabla G (Helmholtz)')

%%

targ.r = [-1:0.001:1; (-1:0.001:1)*0];


[~,~,hess] = helm1d.gshelm(rts,ejs,src,targ);

figure(3)
tiledlayout(1,3);
nexttile
plot(targ.r(1,:),real(hess),targ.r(1,:),imag(hess))
title('\Delta G_S (nonlocal Helmholtz)')

[~,~,hess] = helm1d.gphihelm(rts,ejs,src,targ);

nexttile
plot(targ.r(1,:),real(hess),targ.r(1,:),imag(hess))
title('\Delta G_\phi (or S[G_S])')

k = rts(abs(angle(rts)) == min(abs(angle(rts))));
[~,~,hess] = helm1d.green(k,src,targ);
hess = 1i*hess(:,:,1)./(2*k);

nexttile
plot(targ.r(1,:),real(hess),targ.r(1,:),imag(hess))
title('\Delta G (Helmholtz)')
