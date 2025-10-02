addpath(genpath('~Documents/GitHub/surface-waves'))
addpath(genpath('~/Documents/GitHub/chunkie'))
addpath(genpath('~/Documents/GitHub/fmm3dbie_log_quad'))

% Parameters
 
alpha_int = 4;
beta_int = 1;
[rts_int,ejs_int] = surfwave.capillary.find_roots(alpha_int,beta_int,1);

alpha_ext = 10;
beta_ext = 1;
[rts_ext,ejs_ext] = surfwave.capillary.find_roots(alpha_ext,beta_ext,1);

zpars = [rts_ext;ejs_ext];

zk_int = rts_int(1)
zk_ext = rts_ext(1)

% Creating geometry

rad = 30;
S = geometries.disk([rad,rad],[],[4 3 6],12);
nch = 24;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = chnkr.move([0;0],[0;0],0,rad);
chnkr = sort(chnkr);

figure(1)
plot(S,rand(S.npatches,1))
hold on
plot(chnkr,'x-')


%% Plotting the source

% rhs = -(-alpha_int*zk_ext^3 + beta_int*zk_ext + 1)*exp(1i*zk_ext*S.r(1,:)) ;
% figure(2)
% 
% errs = surf_fun_error(S, rhs); plot(S, log(errs));
% colorbar
% 
% figure(3)
% plot(S,real(rhs)); shading interp
% colorbar


%% Forming matrix for G_S (volume -> volume)

eps = 1e-8;

tic
Gs = surfwave.capillary.gsmatgen(S, zpars, eps);
toc

%% Forming matrix for G_phi (volume -> volume)

tic
Gphi = surfwave.capillary.gphimatgen(S, zpars, eps);
toc

%% Forming matrix for S and D (boundary -> volume)

f3 = chnkr.r(1,:).' / 30;
f4 = zeros(chnkr.npt,1);

fkern1 = @(s,t) surfwave.kern(rts, ejs, s, t,'gs_s');
fkern2 = @(s,t) surfwave.kern(rts, ejs, s, t,'gphi_s');
fkern3 = @(s,t) surfwave.kern(rts, ejs, s, t,'gs_d');
fkern4 = @(s,t) surfwave.kern(rts, ejs, s, t,'gphi_d');

rho1 = - alpha_ext/2*f4;
rho2 = alpha_ext/2*f3;

targobj = [];
targobj.r = S.r(1:2,:);

Ssrho1 = chunkerkerneval(chnkr,fkern1,rho1,targobj);
Sphirho1 = chunkerkerneval(chnkr,fkern2,rho1,targobj);
Dsrho2 = chunkerkerneval(chnkr,fkern3,rho2,targobj);
Dphirho2 = chunkerkerneval(chnkr,fkern4,rho2,targobj);

%% Forming matrices needed for evaluation

[tx,ty] = meshgrid(-60:0.75:60);
targs = [tx(:) ty(:)].';

out = ~chunkerinterior(chnkr,targs);

targs = [targs; 0*tx(:).'] ;
targinfo = [];
targinfo.r = targs(:,out);
ntarg = sum(out)

tic
evalmats = surfwave.capillary.gsevalmat(S, targinfo,zpars, eps);
toc

Gseval = reshape(evalmats(1,:), [S.npts ntarg]);
Gphieval = reshape(evalmats(2,:), [S.npts ntarg]);

load gong.mat
sound(y)

%%

rhs = -(-alpha_int*zk_ext^3 + beta_int*zk_ext + 1)*exp(1i*zk_ext*S.r(1,:)) ;
lhs = -alpha_int/alpha_ext*eye(S.npts) + 0.5*(beta_int - beta_ext*alpha_int/alpha_ext)*Gs.' + (1 - alpha_int/alpha_ext)*Gphi.';

sol = lhs\rhs.';

figure(5)
tiledlayout(2,2)

nexttile
plot(S,real(exp(1i*zk_ext*S.r(1,:)).')); shading interp
colorbar
title("Incident field")
view(0,90)


nexttile
plot(S,real(sol)); shading interp
colorbar
title("Density")
view(0,90)

phizint = 1/2*Gs.'*sol + zk_ext*exp(1i*zk_ext*S.r(1,:)).';

nexttile
plot(S,real(phizint)); shading interp
colorbar
title("\phi_z")
view(0,90)

phiint = Gphi.'*sol + exp(1i*zk_ext*S.r(1,:)).';

nexttile
plot(S,real(phiint)); shading interp
colorbar
title("\phi")
view(0,90)

%%

phizext = tx(:)*0;
xs = tx(:);
phizext(out) = 1/2*Gseval.'*sol + zk_ext*exp(1i*zk_ext*xs(out));
phizext = reshape(phizext, size(tx));

phiext = tx(:)*0;
phiext(out) = Gphieval.'*sol + exp(1i*zk_ext*xs(out));
phiext = reshape(phiext, size(tx));


% s = pcolor(tx,ty,real(uext));
% s.EdgeColor = 'none';
% colorbar

% uinc = tx(:)*0;
% uinc(out) = exp(1i*zk_ext*tx(out)).';
% uinc = reshape(uinc, size(tx));

% s = pcolor(tx,ty,real(uinc));
% s.EdgeColor = 'none';
% colorbar

% phiz = zk_ext*uinc+phizext;
% utotext(utotext == 0) = NaN;

% s = pcolor(tx,ty,real(utotext));
% s.EdgeColor = 'none'; s.FaceColor = 'interp';
% colorbar

figure(8)
tiledlayout(1,2)

nexttile
surf(tx,ty,tx*0,real(phizext),'EdgeColor','none')
hold on
plot(chnkr,'k-','LineWidth',4)
hold on
view(0,90)
axis equal
axis square
plot(S,real(phizint)); shading interp
colorbar
title("\phi_z")
view(0,90)

nexttile
surf(tx,ty,tx*0,real(phiext),'EdgeColor','none')
hold on
plot(chnkr,'k-','LineWidth',4)
hold on
view(0,90)
axis equal
axis square
plot(S,real(phiint)); shading interp
colorbar
title("\phi")
view(0,90)

return

%% check error in the exterior

center = [30;45];
h = 0.01;
xs = (-4*h:h:4*h) + center(1);
ys = (-4*h:h:4*h) + center(2);

targs = [xs.' 0*xs.'+center(2) 0*xs.'; 
         0*ys.'+center(1) ys.' 0*ys.'].';

targinfo = [];
targinfo.r = targs;

tic
evalmats2 = surfwave.capillary.gsevalmat(S, targinfo, zpars, eps);
toc

Gseval2 = reshape(evalmats2(1,:), [S.npts 18]);
Gphieval2 = reshape(evalmats2(2,:), [S.npts 18]);

utest = Gseval2.'*sol;
Sutest = Gphieval2.'*sol;

d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560]/h^2;

lapu = d2*utest(1:9) + d2*utest(10:18);

err2 = 1/2*(alpha_ext*lapu +  beta_ext*utest(5)) + Sutest(5);
err2 = abs(err2)/ abs(max([1/2*alpha_ext*lapu, 1/2*beta_ext*utest(5), Sutest(5)]))


%% check error in the interior

wts = cat(1,S.weights{:});

phitest = real(zk_int*exp(1i*zk_int*S.r(1,:))) ;
Sphitest = real(exp(1i*zk_int*S.r(1,:))) ;
lapphitest = get_surface_laplacian(S,phitest);

err2 = (alpha_int*lapphitest +  beta_int*phitest) + Sphitest;
err2 = abs(err2.*sqrt(wts)) / max([alpha_int*lapphitest, beta_int*phitest, Sphitest]);
err2 = reshape(abs(err2),[],S.npatches);
err2 = log10(max(err2)).';

figure(9)
tiledlayout(2,2)

nexttile
plot(S,real(phizint))
colorbar
view(0,90)
title('\phi_z')

nexttile
errs = surf_fun_error(S, uinttest); plot(S, log10(errs));
colorbar
title('approximation error')
view(0,90)

nexttile
plot(S,err2)
colorbar
view(0,90)
title('control error (plane wave)')

uinttest = real(phizint) ;
Suinttest = real(phiint) ;
lapuinttest = get_surface_laplacian(S,uinttest) ;

err2 = (alpha_int*lapuinttest +  beta_int*uinttest) + Suinttest ;
err2 = abs(err2.*sqrt(wts)) / max([alpha_int*lapuinttest; beta_int*uinttest; Suinttest]);
err2 = reshape(abs(err2),[],S.npatches);
err2 = log10(max(err2)).';

nexttile
plot(S,err2)
colorbar
view(0,90)
title('consistency with PDE')



%% check error against adaptive integration of G(x,y) f(y) 

x = S.r(:,200);
testintgs = integral2(@(y1,y2) integs(a,rts,ejs,x,y1,y2,c1,c2),-30,30,-30,30,'AbsTol',1e-10,'RelTol',1e-10);

err1 = abs(testintgs - Gsf(200)) / abs(Gsf(200))

x = S.r(:,600);
testintgs = integral2(@(y1,y2) integs(a,rts,ejs,x,y1,y2,c1,c2),-30,30,-30,30,'AbsTol',1e-10,'RelTol',1e-10);

err2 = abs(testintgs - Gsf(600)) / abs(Gsf(600))

x = S.r(:,200);
testintgphi = integral2(@(y1,y2) integphi(a,rts,ejs,x,y1,y2,c1,c2),-30,30,-30,30,'AbsTol',1e-10,'RelTol',1e-10);

err3 = abs(testintgphi - Gphif(200) ) / abs(Gphif(200))

x = S.r(:,600);
testintgphi = integral2(@(y1,y2) integphi(a,rts,ejs,x,y1,y2,c1,c2),-30,30,-30,30,'AbsTol',1e-10,'RelTol',1e-10);

err4 = abs(testintgphi - Gphif(600) ) / abs(Gphif(600))