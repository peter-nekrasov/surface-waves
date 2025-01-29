clear 
addpath(genpath('..')) 

zk = 1;
nu = 1/3;

% thetas = 0:pi/6:2*pi-pi/6;
% [targets, ~, ~] = droplet(thetas);
% targets = targets*1.5;

[X,Y] = meshgrid(-6:0.1:6);
targets = [X(:) Y(:)].';
[~,na] = size(targets);

centre = [0. ; 0.];
    
cparams = [];
cparams.maxchunklen = 2;

chnkr = chunkerfunc(@(t) circle(t),cparams);
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
    
coefs = [nu; 0];
opts = [];
opts.sing = 'log';

kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.quad = 'native';
opts2.sing = 'smooth';

fkern1 =  @(s,t) flex2d.kern(zk, s, t, 'free plate first part', coefs);        % build the desired kernel
fkern1bh =  @(s,t) flex2d.kern(zk, s, t, 'free plate first part bh', coefs);        % build the desired kernel
fkern2bh =  @(s,t) flex2d.kern(zk, s, t, 'free plate hilbert bh', coefs);        % build the desired kernel
fkern2 =  @(s,t) flex2d.kern(zk, s, t, 'free plate hilbert', coefs);        % build the desired kernel
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

[~, ~, hess, third, ~] = flex2d.hkdiffgreen(zk, centre, chnkr.r);

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds); % normalization
tauy = (dy./ds);

firstbc = 1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))+...
           coefs(1)/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy));

secondbc = 1./(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny))+...
        (2-coefs(1))/(2*zk^2).*(third(:, :, 1).*(taux.*taux.*nx) + third(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        third(:, :, 3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + third(:, :, 4).*(tauy.*tauy.*ny))+...
        (1-coefs(1)).*(kappa).*(1/(2*zk^2).*(hess(:, :, 1).*taux.*taux + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*tauy.*tauy)-...
        (1/(2*zk^2).*(hess(:, :, 1).*nx.*nx + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*ny.*ny)));

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

ikern1 = @(s,t) flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) flex2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

coupled = chunkerkerneval(chnkr, ikern3, H*rho1, targets(:,out));


start1 = tic;
Dsol = chunkerkerneval(chnkr, ikern1,rho1, targets(:,out)) + coupled +chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

[val,~] = flex2d.hkdiffgreen(zk,centre,targets(:,out));     

utarg(out) = Dsol;
true_sol(out) = 1/(2*zk^2).*val;


%%
utarg = reshape(utarg,size(X));
true_sol = reshape(true_sol,size(X));
uerr = utarg - true_sol;

figure(1);
tiledlayout(1,2)
nexttile
h = pcolor(X,Y,real(utarg));
colorbar
hold on
set(h,'EdgeColor','None'); hold on;

nexttile 
h = pcolor(X,Y,log10(abs(true_sol-utarg)));
colorbar
hold on
set(h,'EdgeColor','None'); hold on;

title("Absolute error (free plate kernel on a starfish)", 'FontSize',16)

