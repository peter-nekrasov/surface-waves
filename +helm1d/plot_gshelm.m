beta = 1+0.5i;
gamma = 2;
[rts,ejs] = helm1d.find_roots(beta,gamma);

src = [];
src.r = [0;0];

targ = [];
targ.r = [-20:0.1:20; (-20:0.1:20)*0];

[val,~] = helm1d.gshelm(rts,ejs,src,targ);

figure(1)
tiledlayout(1,3);
nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G_S (nonlocal Helmholtz)')
ylim([-0.3 0.3])


[val,~] = helm1d.gphihelm(rts,ejs,src,targ);

nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G_\phi (or S[G_S])')
ylim([-0.3 0.3])

k = rts(abs(angle(rts)) == min(abs(angle(rts))));
[val,~] = helm1d.green(k,src,targ);

nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G (Helmholtz)')
ylim([-0.3 0.3])

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
[~,grad] = helm1d.green(k,src,targ);
grad = grad(:,:,1);

nexttile
plot(targ.r(1,:),real(grad),targ.r(1,:),imag(grad))
title('\nabla G (Helmholtz)')

%%

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
hess = hess(:,:,1);

nexttile
plot(targ.r(1,:),real(hess),targ.r(1,:),imag(hess))
title('\Delta G (Helmholtz)')
