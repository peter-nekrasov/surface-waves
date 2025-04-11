addpath(genpath('~Documents/GitHub/surface-waves'))
addpath(genpath('~Documents/GitHub/chunkie'))

zk = 8;
cparams = [];
cparams.maxchunklen = 4/zk; 

chnkr = chunkerfunc(@(t) ellipse(t,2,1), cparams);
chnkr = chnkr.sort();

centre = [0.1; 0.2];

alpha = 1;
beta = -2;
gamma = 10;

[rts,ejs] = surfwave.find_roots(alpha,beta,gamma);

angles = angle(rts);

% Find the index of the smallest angle
[~, idx] = min(abs(angles));

% Get the complex number with the smallest angle
k_real = rts(idx)

% PLOTTING GEOMETRY

figure(1)                                                   % plot the chunker-object (supposed to be a starfish3arms centered at 1 with radius 1)
clf
axis equal
% polarplot(ones(1,200), 'k' )
% set(gca, 'ThetaTick', (0:pi/2:3*pi/2)*180/(pi))
% set(gca, 'ThetaTickLabel', {'0','\pi/2','\pi','3\pi/2'})
plot(chnkr, '-kx','LineWidth',1)
hold on
xlim([-3 3])
ylim([-3 3])
title('Scattering object')

% HELMHOLTZ NEUMANN

fkern = @(s,t) surfwave.kern(rts,ejs, s, t,'gs_d');
fkern3 = @(s,t) surfwave.kern(rts, ejs, s, t,'gphi_d');

sysmat = chunkermat(chnkr,fkern);

lhs = eye(chnkr.k*chnkr.nch) + sysmat;

rhs = chnkr.r(1,:).';
sol = lhs\rhs;

[tx,ty] = meshgrid(-3:0.06:3);
targs = [tx(:) ty(:)].';

in = chunkerinterior(chnkr,targs);
out = ~in;

usol = out*0;

usol(out) = chunkerkerneval(chnkr,fkern,sol,targs(:,out)); 
usol(~out) = NaN;
usol = reshape(usol,size(tx));

figure(2)
t = pcolor(tx,ty,real(usol));
t.EdgeColor = 'None';
hold on
colorbar
plot(chnkr, '-k','LineWidth',2)
title('\Re(\sigma)')

usol = out*0;

usol(out) = chunkerkerneval(chnkr,fkern3,sol,targs(:,out)); 
usol(~out) = NaN;
usol = reshape(usol,size(tx));

figure(3)
t = pcolor(tx,ty,real(usol));
t.EdgeColor = 'None';
hold on
colorbar
plot(chnkr, '-k','LineWidth',2)
title('\Re(S[\sigma])')


%% error in exterior PDE

h = 0.05;
center = 1.5;
xs = center-4*h:h:center+4*h;
ys = center-4*h:h:center+4*h;
hortargs = [xs; xs*0+center];
vertargs = [ys*0+center; ys];

horvals = chunkerkerneval(chnkr,fkern,sol,hortargs); 
vertvals = chunkerkerneval(chnkr,fkern,sol,vertargs); 

sus = chunkerkerneval(chnkr,fkern3,sol,[center;center]); 

us = zeros(9);
us(5,:) = horvals;
us(:,5) = vertvals;

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

f1 = 0.5*alpha*sum(lap.*us,'all') + 0.5*beta*us(5,5) + gamma*sus;
err = abs(f1) / max([0.5*alpha*sum(lap.*us,'all'), 0.5*beta*us(5,5), gamma*sus]) 

%% error in Dirichlet BC

h = 0.25;

bdy_r = chnkr.r(:,4,4);
bdy_n = chnkr.n(:,4,4);

theta = atan2(bdy_n(2),bdy_n(1));

r1 = bdy_r(1) + h*cos(theta);
r2 = bdy_r(2) + h*sin(theta);

targ = [r1;r2];

us = chunkerkerneval(chnkr,fkern,sol,targ); 

err = abs(us - bdy_r(1)) / abs(bdy_r(1))
