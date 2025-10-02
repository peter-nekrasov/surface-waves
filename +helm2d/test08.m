%% Finite difference test for derivatives of axisymmetric Greens functions

addpath(genpath("~/Documents/GitHub/chunkie"))

h = 0.05;
xs = 2:h:(2+8*h);
alpha = 0.5;
beta = 1.5;
gamma = 2.5;

n = 0;

src = [];
src.r = [1;0];
targ = [];
targ.r = [xs; xs*0];

[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);

[gs,gsgrad,gshess] = helm2d.gsaxisym(rts,ejs,n,src,targ);

[gphi,gphigrad,gphihess] = helm2d.gphiaxisym(rts,ejs,n,src,targ);

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

err1 = d2*gs + 1/xs(5)*d1*gs - n^2/xs(5)^2*gs(5) - gshess(5)
err2 = d1*gs - gsgrad(5)
err3 = d2*gphi + 1/xs(5)*d1*gphi - n^2/xs(5)^2*gphi(5) - gphihess(5)
err4 = d1*gphi - gphigrad(5)

% passed