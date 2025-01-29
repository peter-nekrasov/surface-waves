clear 


% zk = 1;
% 
% %
% 
% zk = 0.1;               % our k (wave number)

zk = 18;
nu = 1/3;
cparams = [];

cparams.eps = 1e-6;
cparams.nover = 0;
cparams.maxchunklen = 0.2;       % setting a chunk length helps when the
                                    % frequency is known'
%cparams.maxchunklen = 0.4;


chnkr = chunkerfunc(@(t) ellipse(t), cparams);

coefs = [nu; 0];
opts = [];
opts.sing = 'log';
fkern =  @(s,t) flex2d.kern(zk, s, t, 'free plate first part', coefs);        % build the desired kernel


fkern1 = @(s,t) flex2d.kern(zk, s, t, 'free plate hilbert subtract', coefs);                   % hilbert subtraction kernels in K11
fkern2 = @(s,t) flex2d.kern(zk, s, t, 'free plate coupled hilbert', coefs);   

hilbert = @(s,t) chnk.lap2d.kern(s, t, 'hilb');
double = @(s,t) chnk.lap2d.kern(s,t, 'd');

fkern3 = @(s,t) flex2d.kern(zk, s, t, 'free plate K21 first part', coefs);                     % singularity subtration kernel in K21 (including swapping its Asmyptotics expansions)

fkern4 = @(s,t) flex2d.kern(zk, s, t, 'free plate K21 second part', coefs);                    % kernels in K21 needs to multiply by curvature

fkern5 = @(s,t) flex2d.kern(zk, s, t, 'free plate K21 hilbert part', coefs);                   % kernels in K21 coupled with hilbert transforms and needs to multiply by curvature

fkern6 = @(s,t) flex2d.kern(zk, s, t, 'free plate K22 second part', coefs);                    % kernels in K22 needs to multiply by curvature


start = tic;
sysmat = chunkermat(chnkr,fkern, opts);
sysmat1 = chunkermat(chnkr, fkern1, opts);
sysmat2 = chunkermat(chnkr, fkern2, opts);

K21 = chunkermat(chnkr, fkern3, opts);

K21second = chunkermat(chnkr, fkern4, opts);
K21hilbert = chunkermat(chnkr, fkern5, opts);

K22second = chunkermat(chnkr, fkern6, opts);

D = chunkermat(chnkr, double, opts);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

opts2 = [];
opts2.sing = 'pv';

start = tic;
H = chunkermat(chnkr, hilbert, opts2);                                              % Assemble hilbert transforms
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)



kappa = signed_curvature(chnkr);
kappa = kappa(:);


hilb = sysmat1*H - ((1+nu)/2).*(D*D)- ((1+nu)*nu/2).*(D*D);
hilb2 = sysmat2*H ;

mat1 =  sysmat(1:2:end, 1:2:end);
mat4 =  sysmat(2:2:end, 2:2:end);


sysmat(1:2:end, 1:2:end) = mat1 + hilb;
sysmat(2:2:end, 1:2:end) = K21 +  hilb2 + (kappa).*(K21hilbert*H + K21second);
sysmat(2:2:end, 2:2:end) = mat4 + (kappa).*(K22second);


A = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for exterior problem)

M = kron(eye(chnkr.npt), A);

lhs =  M + sysmat ;

zkimag = (1i)*zk;
[~, ~, hess, third, ~] = flex2d.helmdiffgreen(zk, [0;0], chnkr.r);
[~, ~, hessK, thirdK, ~] = flex2d.helmdiffgreen(zkimag, [0;0], chnkr.r);

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

firstbc = 1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))-...
           1/(2*zk^2).*(hessK(:, :, 1).*(nx.*nx) + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*(ny.*ny))+...
           coefs(1)/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))-...
           coefs(1)/(2*zk^2).*(hessK(:, :, 1).*(taux.*taux) + hessK(:, :, 2).*(2*taux.*tauy) + ...
           hessK(:, :, 3).*(tauy.*tauy));

secondbc = 1./(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nx.*nx.*nx) + thirdK(:, :, 2).*(3*nx.*nx.*ny)+...
        thirdK(:, :, 3).*(3*nx.*ny.*ny) + thirdK(:, :, 4).*(ny.*ny.*ny)) +...
        (2-coefs(1))/(2*zk^2).*(third(:, :, 1).*(taux.*taux.*nx) + third(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        third(:, :, 3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + third(:, :, 4).*(tauy.*tauy.*ny)) - ...
        (2-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(taux.*taux.*nx) + thirdK(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        thirdK(:, :, 3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + thirdK(:, :, 4).*(tauy.*tauy.*ny)) +...
        (1-coefs(1)).*(kappa).*(1/(2*zk^2).*(hess(:, :, 1).*taux.*taux + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*tauy.*tauy)-...
         1/(2*zk^2).*(hessK(:, :, 1).*taux.*taux + hessK(:, :, 2).*(2*taux.*tauy) + ...
        hessK(:, :, 3).*tauy.*tauy)-...
        (1/(2*zk^2).*(hess(:, :, 1).*nx.*nx + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*ny.*ny)-...
       1/(2*zk^2).*(hessK(:, :, 1).*nx.*nx + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*ny.*ny)));

[nt, ~] = size(sysmat);

rhs = zeros(nt, 1); 
rhs(1:2:end) = firstbc ; 
rhs(2:2:end) = secondbc;

tic
%sol = gmres(lhs, rhs, [], 1e-13, 200);
sol = lhs\rhs;
toc;

rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);        


xs = -4:0.05:4;                                    % generate some targets
ys = -4:0.05:4;
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

tic
in = chunkerinterior(chnkr, targets); 
out = ~in; 
toc


ikern1 = @(s,t) flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) flex2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

coupled = chunkerkerneval(chnkr, ikern3, H*rho1, targets(:, out));


start1 = tic;
Dsol = chunkerkerneval(chnkr, ikern1,rho1, targets(:, out)) + coupled +chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

true_sol = zeros(na, 1);
utarg = zeros(na, 1);

[val,~] = flex2d.helmdiffgreen(zk,[0;0],targets(:,out));        % Hankel part

zkimag = (1i)*zk;
[valK,~] = flex2d.helmdiffgreen(zkimag, [0;0], targets(:,out));    % modified bessel K part

trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;


utarg(out) = Dsol;
true_sol(out) = trueval;


uerr = utarg - true_sol;
uerr = reshape(uerr,size(X));
figure
h = pcolor(X,Y,log10(abs(uerr)));
set(h,'EdgeColor','None'); hold on;
title("Absolute error (free plate kernel on a ellipse k = 18)", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);

colorbar

