%% Place to check axisymmetric greens stuff

h = 0.05;
xs = 2:h:(2+8*h);
alpha = 0.5;
beta = 1.5;
gamma = 2.5;

n = 0;

src = [];
src.r = [1 2 3;0 0 0];
targ = [];
targ.r = [5 6 7;0 0 0];

[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
[gs,grad,hess] = helm2d.gsaxisym(rts,ejs,n,src,targ);
