close all

alpha = 1;
beta = 0;
gamma = 1;
[rts,ejs] = surfwave.find_roots(alpha,beta,gamma);

rts

[tx, tz] = meshgrid(-20:0.01:20,-5:0.01:0);

val1 = helm1d.domain_gf(rts,ejs,tx-7*pi/2,tz);
val2 = helm1d.domain_gf(rts,ejs,tx+7*pi/2,tz);

% plotting a trapped mode
figure(1)
s = pcolor(tx,tz,real(val2+val1));
s.EdgeColor = 'None';
colorbar

xs = (-20:0.005:20);
val1 = helm1d.domain_gf(rts,ejs,xs-7*pi/2,0);
val2 = helm1d.domain_gf(rts,ejs,xs+7*pi/2,0);

figure(2)
plot(xs,real(val1+val2),xs,imag(val1+val2))

val1 = helm1d.domain_gf(rts,ejs,xs-7*pi/2,-0.1);
val2 = helm1d.domain_gf(rts,ejs,xs+7*pi/2,-0.1);

figure(3)
plot(xs,real(val1+val2),xs,imag(val1+val2))

figure(4)
zs = (-20:0.005:0);
val1 = helm1d.domain_gf(rts,ejs,-7*pi/2,zs);
val2 = helm1d.domain_gf(rts,ejs,7*pi/2,zs);
plot(zs,real(val1+val2),zs,imag(val1+val2))

