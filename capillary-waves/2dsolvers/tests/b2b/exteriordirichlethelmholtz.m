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

fkern = @(s,t) chnk.helm2d.kern(zk, s, t,'d');

sysmat = chunkermat(chnkr,fkern);

lhs = 0.5*eye(chnkr.k*chnkr.nch) + sysmat;

nx = chnkr.n(1,:).' ; ny = chnkr.n(2,:).';

src = [];
src.r = centre;

firstbc = chnk.helm2d.green(zk,src.r,chnkr.r);

rhs = firstbc;
sol = lhs\rhs;

[tx,ty] = meshgrid(-3:0.05:3);
targs = [tx(:) ty(:)].';

in = chunkerinterior(chnkr,targs);
out = ~in;

usol = out*0;
truesol = out*0;

usol(out) = chunkerkerneval(chnkr,fkern,sol,targs(:,out)); 
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

return