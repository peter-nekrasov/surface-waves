%% Place to check axisymmetric greens stuff

h = 0.05;
xs = 2:h:(2+8*h);
alpha = 0.5;
beta = 1.5;
gamma = 2.5;

n = 0;

targ.r = [10.0166;0];
src.r = [10;0];

[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
[gs,grad,hess] = helm2d.gsaxisym(rts,ejs,n,src,targ);
