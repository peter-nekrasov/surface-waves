% CREATING GEOMETRY  - remember to check kappa prime in kern.m

addpath(genpath('~Documents/GitHub/chunkie'))

zk = 8;
cparams = [];
cparams.maxchunklen = 8/zk; 

theta = 0;

chnkr = chunkerfunc(@(t) ellipse(t), cparams);
chnkr = chnkr.sort();

centre = [0.1; 0.2];


% PLOTTING GEOMETRY

figure(1)                                                   % plot the chunker-object (supposed to be a starfish3arms centered at 1 with radius 1)
clf
axis equal
% polarplot(ones(1,200), 'k' )
% set(gca, 'ThetaTick', (0:pi/2:3*pi/2)*180/(pi))
% set(gca, 'ThetaTickLabel', {'0','\pi/2','\pi','3\pi/2'})
plot(chnkr, '-k','LineWidth',2)
hold on
xlim([-3 3])
ylim([-3 3])
title('Scattering object')

% HELMHOLTZ NEUMANN

fkern = @(s,t) chnk.helm2d.kern(zk, s, t,'sprime');
fkern2 = @(s,t) chnk.helm2d.kern(zk, s, t,'s');

sysmat = chunkermat(chnkr,fkern);

lhs = -0.5*eye(chnkr.k*chnkr.nch) + sysmat;

nx = chnkr.n(1,:).' ; ny = chnkr.n(2,:).';

src = [];
src.r = centre;

[~,grad] = chnk.helm2d.green(zk,src.r,chnkr.r);
firstbc = grad(:,1).*nx + grad(:,2).*ny;

rhs = firstbc;
rhs = chnkr.r(1,:).';

sol = lhs\rhs;

[tx,ty] = meshgrid(-3:0.05:3);
targs = [tx(:) ty(:)].';

in = chunkerinterior(chnkr,targs);
out = ~in;

usol = out*0;
truesol = out*0;

usol(out) = chunkerkerneval(chnkr,fkern2,sol,targs(:,out)); 
truesol(out) =  chnk.helm2d.green(zk,src.r,targs(:,out));

usol(~out) = NaN;
truesol(~out) = NaN;

usol = reshape(usol,size(tx));
truesol = reshape(truesol,size(tx));

err = usol - truesol;

figure(2)
t = pcolor(tx,ty,real(usol));
t.EdgeColor = 'None';
hold on
colorbar
plot(chnkr, '-k','LineWidth',2)
title('\Re(u)')

figure(3)
t = pcolor(tx,ty,real(truesol));
t.EdgeColor = 'None';
hold on
colorbar
plot(chnkr, '-k','LineWidth',2)
title('\Re(u_{ref})')


figure(4)
t = pcolor(tx,ty,log10(abs(err)));
t.EdgeColor = 'None';
hold on
colorbar
plot(chnkr, '-k','LineWidth',2)
title('Log_{10} error')


%% 

h = 0.00001;

bdy_r = chnkr.r(:,1,1);
bdy_n = chnkr.n(:,1,1);
bdy_tau = chnkr.d(:,1,1);

theta = atan2(bdy_n(2),bdy_n(1));

d1 = [-49/20	6	-15/2	20/3	-15/4	6/5	-1/6]/h;

pts = h:h:7*h;

r1 = bdy_r(1) + pts*cos(theta);
r2 = bdy_r(2) + pts*sin(theta);

targs = [r1;r2];

figure(1)
hold on
plot(r1,r2,'x')
xlim([min(r1)-2*h max(r1)+2*h])
ylim([min(r2)-2*h max(r2)+2*h])
axis equal

us = chunkerkerneval(chnkr,fkern2,sol,targs); 

dudn = d1*us

err = abs(bdy_r(1) - dudn) / abs(dudn)



return