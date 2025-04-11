%% Using trapezoid rule to integrate the axisymmetric G against some rhs f
% adaptive integration too slow :(

h = 0.01;
rs = 8:h:(8+8*h);
center = 8+4*h;
targ = [];
targ.r = [rs; 0*rs];

alpha = 0.5;
beta = 1.5;
gamma = 2.5;

n = 2;

[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);

f0 = f(center);

src = []; src.r = [0:h:20; 0*(0:h:20)];

gsmat = helm2d.gsaxisym(rts,ejs,n,src,targ);

rhs = f(0:h:20).';
usol = -gsmat*rhs.*h;

targ = []; targ.r = [center; 0];
gphimat = helm2d.gphiaxisym(rts,ejs,n,src,targ);
Susol = -gphimat*rhs.*h;

%%

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

f0
f1 = alpha/2*d2*usol + alpha/2/rs(5)*d1*usol - alpha/2*n^2/rs(5)^2*usol(5) + beta/2*usol(5) + gamma*Susol 

err = f0 - f1 %  

function val = integrand(rts,ejs,n,targ,rhop)
    src = [];
    src.r = [rhop; 0];
    gs = helm2d.gsaxisym(rts,ejs,n,src,targ);
    val = -gs.*f(rhop);
end

function val = integrand2(rts,ejs,targ,rhop)
    src = [];
    src.r = [rhop; 0];
    gphi = helm2d.gphihelm(rts,ejs,n,src,targ);
    val = -gphi.*f(rhop);
end

function val = f(r)
    val = exp(-((r-10).^2)/2);
end

