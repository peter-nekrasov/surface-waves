clear
close all 
addpath(genpath('../..'))

[rs,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);
rs = rs/8 + 1/8;

zk = 2;
R = sqrt(2);

[X1,X2]=meshgrid(-8:0.1:8);
u = 0*X1;
true_sol = 0*X1;
RX = sqrt(X1.^2+X2.^2);
out = (RX > 1.1*R);
targs = [X1(out) X2(out)].';

% calculate fourier modes 

N = 100;
dt = 2*pi/N;
thetas = dt:dt:2*pi;
rho1 = thetas*0;
rho2 = thetas*0;

src = [0.3;-0.3];
targ = [R;0];

for n = -N/4:N/4

    k11s = [];
    k12s = [];
    k21s = [];
    k22s = [];
    
    for i = 1:16
        R2 = R+rs(i);
        k11s = [k11s; integral(@(t) exp(1i*n*t).*K11(zk,R,R2,t),0,2*pi)];
        k12s = [k12s; integral(@(t) exp(1i*n*t).*K12(zk,R,R2,t),0,2*pi)];
        k21s = [k21s; integral(@(t) exp(1i*n*t).*K21(zk,R,R2,t),0,2*pi)];
        k22s = [k22s; integral(@(t) exp(1i*n*t).*K22(zk,R,R2,t),0,2*pi)];
    end
    
    k11 = pols.'*v2c*k11s;
    k12 = pols.'*v2c*k12s;
    k21 = pols.'*v2c*k21s;
    k22 = pols.'*v2c*k22s;
    
    lhs = [k11, k12; k21, k22];
    
    bc1 = integral(@(t) exp(1i*n*t).*BC1(zk,sqrt(src(1)^2+src(2)^2),R,t),0,2*pi);
    bc2 = integral(@(t) exp(1i*n*t).*BC2(zk,sqrt(src(1)^2+src(2)^2),R,t),0,2*pi);
    rhs = exp(-1i*n*atan2(src(2),src(1)))*[bc1;bc2];
    
    sol = lhs \ rhs;
    
    rho1 = rho1+sol(1)*exp(1i*n*thetas);
    rho2 = rho2+sol(2)*exp(1i*n*thetas);

end

srcs = R*[cos(thetas); sin(thetas)];

[~,~,hess,third] = flex2d.hkdiffgreen(zk,srcs,targs);

[~,nt] = size(targs);

nx = cos(thetas).';
ny = sin(thetas).';

taux = -sin(thetas).';
tauy = cos(thetas).';

nx = repmat(nx.',[nt 1]);
ny = repmat(ny.',[nt 1]);
taux = repmat(taux.',[nt 1]);
tauy = repmat(tauy.',[nt 1]);

K1 = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
   third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny))) - ...
   (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
   third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

K2 = -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
      (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

u(out) = K1*rho1.'*dt/(2*pi) + K2*rho2.'*dt/(2*pi);
 
[G0,~] = flex2d.hkdiffgreen(zk,src,targs);
true_sol(out) = G0/(2*zk^2);

%%


set(groot, 'DefaultSurfaceEdgeColor', 'none')
tiledlayout(1,3);
nexttile
pcolor(X1,X2,real(u));
hold on
pos = [-R -R 2*R 2*R];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor','white')
title('BIE/Fourier Solver')
colorbar

nexttile
pcolor(X1,X2,real(true_sol));
rectangle('Position',pos,'Curvature',[1 1],'FaceColor','white')
title('Analytic Solution')
colorbar

nexttile
err = abs(u - true_sol) ;
pcolor(X1,X2,log10(abs(err)));
rectangle('Position',pos,'Curvature',[1 1],'FaceColor','white')
title('Absolute Error')
colorbar

function val = K11(zk,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,~,third] = flex2d.hkdiffgreen(zk,src,targ);

    nx = cos(t);
    ny = sin(t);

    taux = -sin(t);
    tauy = cos(t);
    
    val = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
         third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny))) - ...
         (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
          third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

end

function val = K12(zk,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,hess] = flex2d.hkdiffgreen(zk,src,targ);

    nx = cos(t);
    ny = sin(t);

    taux = -sin(t);
    tauy = cos(t);
    
    val = -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

end

function val = K21(zk,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,~,~,forth] = flex2d.hkdiffgreen(zk,src,targ);

    nx = cos(t);
    ny = sin(t);

    taux = -sin(t);
    tauy = cos(t);

    nxtarg = 1;
    nytarg = 0;
    
    val = -(1/(2*zk^2).*(forth(:, :, 1).*(nx.*nx.*nx.*nxtarg) + forth(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
      forth(:, :, 3).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + forth(:, :, 4).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
      forth(:, :, 5).*(ny.*ny.*ny.*nytarg)) ) - ...
      (3/(2*zk^2).*(forth(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ forth(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
      forth(:, :, 3).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
      forth(:, :, 4).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) +...
      forth(:, :, 5).*(ny.*tauy.*tauy.*nytarg)));

end

function val = K22(zk,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,~,third] = flex2d.hkdiffgreen(zk,src,targ);

    nx = cos(t);
    ny = sin(t);

    taux = -sin(t);
    tauy = cos(t);

    nxtarg = 1;
    nytarg = 0;
    
    val = -(1/(2*zk^2).*(third(:,:, 1).*(nx.*nx.*nxtarg) +third(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + third(:, :, 3).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         third(:, :,4).*(ny.*ny.*nytarg))) + ...
         (1/(2*zk^2).*(third(:,:, 1).*(taux.*taux.*nxtarg) +third(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + third(:, :, 3).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         third(:, :,4).*(tauy.*tauy.*nytarg)));

end


function val = BC1(zk,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    val = flex2d.hkdiffgreen(zk,src,targ);
    
    val = 1/(2*zk^2).*val;

end

function val = BC2(zk,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,grad] = flex2d.hkdiffgreen(zk,src,targ);

    nx = 1;
    ny = 0;
    
    val = 1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny );

end