beta = 2+0.4i;
gamma = 1;
[rts,ejs] = helm1d.find_roots(beta,gamma);

src = [];
src.r = [50;0];

targ = [];
targ.r = [0:0.01:100; (0:0.01:100)*0];

[val,~] = helm1d.gshelm(rts,ejs,src,targ);

tiledlayout(1,2);
nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G_S (nonlocal Helmholtz)')
ylim([-0.3 0.3])

k = rts(1);
[val,~] = helm1d.green(k,src,targ);
val = - val;

nexttile
plot(targ.r(1,:),real(val),targ.r(1,:),imag(val))
title('G (Helmholtz)')
ylim([-0.3 0.3])