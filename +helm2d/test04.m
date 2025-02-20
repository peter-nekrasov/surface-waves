%% Finite difference test for axisymmetric Green's functions
% see notes for the radial ODE that is being tested

h = 0.05;
xs = 2:h:(2+8*h);
alpha = 0.5;
beta = 1.5;
gamma = 2.5;

n = 2;

src = [];
src.r = [1;0];
targ = [];
targ.r = [xs; xs*0];

[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);

gs = helm2d.gsaxisym(rts,ejs,n,src,targ);

targ.r = [xs(5); 0];
gphi = helm2d.gphiaxisym(rts,ejs,n,src,targ);

% finite difference stencils
d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;

d2 = zeros(1, 9);
d2(1) = -1/560;
d2(2) = 8/315;
d2(3) = -1/5;
d2(4) = 8/5;
d2(5) = -205/72;
d2(6) = 8/5;
d2(7) = -1/5;
d2(8) = 8/315;
d2(9) = -1/560;
d2 = d2 / h^2;

err = abs( alpha/2*d2*gs + alpha/2/xs(5)*d1*gs - alpha/2*n^2/xs(5)^2*gs(5) + beta/2*gs(5) + gamma*gphi ) ./ max(abs(gs(:))) 

% passed