close all

alpha = 1;
beta = 600+0.2i;
gamma = 1;
[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
k = rts(abs(angle(rts)) == min(abs(angle(rts))));

src = [];
src.r = [0;0];

targ = [];
targ.r = [-1:0.01:1; (-1:0.01:1)*0];

[val,~] = helm2d.gshelm(rts,ejs,src,targ);

figure(1)
tiledlayout(1,3);
nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G_S (nonlocal Helmholtz)')


[val] = helm2d.gphihelm(rts,ejs,src,targ);

nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G_\phi (or S[G_S])')

[val,~] = chnk.helm2d.green(k,src.r,targ.r);

nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G (Helmholtz)')

%%

[~,grad] = helm2d.gshelm(rts,ejs,src,targ);
grad = grad(:,:,1);

figure(2)
tiledlayout(1,3);
nexttile
plot(targ.r(1,:),real(grad),targ.r(1,:),imag(grad))
title('\nabla G_S (nonlocal Helmholtz)')


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
