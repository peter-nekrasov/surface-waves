clear 
addpath(genpath('..')) 

beta = 3;
gamma = -1;
nu = 1/3;

[rts,ejs] = flex2d.find_roots(beta,gamma);

[X,Y] = meshgrid(-6:0.1:6);
targets = [X(:) Y(:)].';
[~,na] = size(targets);

% thetas = 0:pi/6:2*pi-pi/6;
% [targets, ~, ~] = droplet(thetas);
% targets = targets*1.5;
centre = [0.8 ; 0.5];
    
cparams = [];
cparams.maxchunklen = 2;

chnkr = chunkerfunc(@(t) droplet(t),cparams);
chnkr = chnkr.sort();

figure(1)                                                   % plot the chunker-object (supposed to be a circle centered at 1 with radius 1)
clf
scatter(targets(1,:),targets(2,:),36,'filled')
hold on 
scatter(centre(1),centre(2),36,'filled')
hold on
plot(chnkr, '-x')
hold on
quiver(chnkr)
hold on 
axis equal
drawnow
    
coefs = {nu, rts, ejs};

opts = [];
opts.sing = 'log';

kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.quad = 'native';
opts2.sing = 'smooth';

fkern1 =  @(s,t) flex2d.kern(1, s, t, 'free plate first part fg', coefs);        % build the desired kernel
fkern1bh =  @(s,t) flex2d.kern(1, s, t, 'free plate first part bh', coefs);        % build the desired kernel
fkern2bh =  @(s,t) flex2d.kern(1, s, t, 'free plate hilbert bh', coefs);        % build the desired kernel
fkern2 =  @(s,t) flex2d.kern(1, s, t, 'free plate hilbert fg', coefs);        % build the desired kernel
double = @(s,t) lap2d.kern(s,t,'d',coefs);
hilbert = @(s,t) lap2d.kern(s,t,'hilb',coefs);

sysmat1 = chunkermat(chnkr,fkern1, opts);
sysmat1bh = chunkermat(chnkr,fkern1bh, opts2);
sysmat2bh = chunkermat(chnkr,fkern2bh, opts2);
sysmat2 = chunkermat(chnkr,fkern2, opts);

D = chunkermat(chnkr, double, opts);

opts3 = [];
opts3.sing = 'pv';

H = chunkermat(chnkr, hilbert, opts3);     

% Perform diagonal replacement for smooth quads here

sysmat1bh(isnan(sysmat1bh)) = 0;
sysmat1bh(2:2:end,1:2:end) = sysmat1bh(2:2:end,1:2:end) + diag((-3+3*nu)/(8*pi)*kappa.^2.*chnkr.wts(:));
sysmat1 = sysmat1 + sysmat1bh;

sysmat2bh(isnan(sysmat2bh)) = 0;
sysmat2 = sysmat2 + sysmat2bh;

sysmat2(1:2:end,1:2:end) = sysmat2(1:2:end,1:2:end)*H  - 2*((1+nu)/2)^2*D*D;
sysmat2(2:2:end,1:2:end) = sysmat2(2:2:end,1:2:end)*H;

sysmat = sysmat1 + sysmat2 ;

D = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for exterior problem)
D = kron(eye(chnkr.npt), D);

lhs =  D + sysmat;

[~, ~, hess, third, ~] = flex2d.fggreen(centre, chnkr.r, rts, ejs);

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds); % normalization
tauy = (dy./ds);

firstbc = (hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))+...
           nu.*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy));

secondbc = (third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny))+...
        (2-nu).*(third(:, :, 1).*(taux.*taux.*nx) + third(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        third(:, :, 3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + third(:, :, 4).*(tauy.*tauy.*ny))+...
        (1-nu).*(kappa).*((hess(:, :, 1).*taux.*taux + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*tauy.*tauy)-...
        ((hess(:, :, 1).*nx.*nx + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*ny.*ny)));

[nt, ~] = size(sysmat);

rhs = zeros(nt, 1); 
rhs(1:2:end) = firstbc ; 
rhs(2:2:end) = secondbc;

sol = lhs\rhs;

rho1 = sol(1:2:end);                              
rho2 = sol(2:2:end);       

in = chunkerinterior(chnkr, targets); 
out = ~in; 

true_sol = zeros(na, 1);
utarg = zeros(na, 1);

ikern1 = @(s,t) flex2d.kern(1, s, t, 'free plate eval first fg', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) flex2d.kern(1, s, t, 'free plate eval second fg',coefs);
ikern3 = @(s,t) flex2d.kern(1, s, t, 'free plate eval first hilbert fg',coefs);

coupled = chunkerkerneval(chnkr, ikern3, H*rho1, targets(:, out));

start1 = tic;
Dsol = chunkerkerneval(chnkr, ikern1,rho1,  targets(:, out)) + coupled +chunkerkerneval(chnkr, ikern2, rho2,  targets(:, out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

[val,~] = flex2d.fggreen(centre,targets(:,out),rts,ejs);     

utarg(out) = Dsol;
true_sol(out) = val;

%%

uerr = utarg - true_sol;
uerr = reshape(utarg,size(X));
figure
h = pcolor(X,Y,real(uerr));
colorbar
hold on
set(h,'EdgeColor','None'); hold on;
title("Absolute error (free plate kernel on a starfish)", 'FontSize',16)
plot(chnkr,'-x');

rmpath(genpath('..')) 
