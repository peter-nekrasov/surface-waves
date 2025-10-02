%% Using adaptive integration to integrate G against some rhs f

h = 0.05;
xs = 1:h:(1+8*h);
ys = 1:h:(1+8*h);
center = 1+4*h;
hortargs = [xs; xs*0+center];
vertargs = [ys*0+center; ys];

alpha = 0.5;
beta = -1.5;
gamma = 2.5;

[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);

f0 = f(center,center);

horgs = zeros(1,9);
vergs = zeros(9,1);

for ii = 1:9

    x = hortargs(:,ii);
    horgs(ii) = integral2(@(y1,y2) integrand(rts,ejs,x,y1,y2),-20,20,-20,20);

    x = vertargs(:,ii);
    vergs(ii) = integral2(@(y1,y2) integrand(rts,ejs,x,y1,y2),-20,20,-20,20);

    disp(ii)
end

x = [center;center];
gphi = integral2(@(y1,y2) integrand2(rts,ejs,x,y1,y2),-20,20,-20,20);
 

%%

gs = zeros(9);
gs(5,:) = horgs;
gs(:,5) = vergs;

% finite difference stencils
d2 = zeros(9, 1);
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

lap = zeros(9,9);
lap(:,5) = d2;
lap(5,:) = lap(5,:) + d2.';

f0
f1 = 0.5*alpha*sum(lap.*gs,'all') + 0.5*beta*gs(5,5) + gamma*gphi

err = f0 - f1 % passed 

function val = integrand(rts,ejs,x,y1,y2)
    src = [];
    src.r = [y1(:) y2(:)].';
    targ = [];
    targ.r = x;
    gs = helm2d.gshelm(rts,ejs,src,targ);
    gs = reshape(gs,size(y1));
    val = -gs.*f(y1,y2);
end

function val = integrand2(rts,ejs,x,y1,y2)
    src = [];
    src.r = [y1(:) y2(:)].';
    targ = [];
    targ.r = x;
    gphi = helm2d.gphihelm(rts,ejs,src,targ);
    gphi = reshape(gphi,size(y1));
    val = -gphi.*f(y1,y2);
end

function val = f(y1,y2)
    val = exp(-(y1.^2+y2.^2)/8);
end

