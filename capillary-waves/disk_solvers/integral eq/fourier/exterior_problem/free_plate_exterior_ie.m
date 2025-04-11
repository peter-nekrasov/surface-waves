clear
close all 
addpath(genpath('../..'))

[rs,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);
rs = rs/8 + 1/8;

R = sqrt(2);
zk = 2;
nu = 0.3;
beta = (1+nu)/2;

[X1,X2]=meshgrid(-8:0.1:8);
u = 0*X1;
true_sol = 0*X1;
RX = sqrt(X1.^2+X2.^2);
out = (RX > 1.05*R);
targs = [X1(out) X2(out)].';

% calculate fourier modes 

N = 100;
dt = 2*pi/N;
thetas = dt:dt:2*pi;
rho1 = thetas*0;
Hrho1 = thetas*0;
rho2 = thetas*0;

src = [0.3;-0.3];
targ = [R;0];

for n = -N/2:N/2

    k11s = [];
    k12s = [];
    k21s = [];
    k22s = [];
    
    for i = 1:16
        R2 = R+rs(i);
        k11s = [k11s; integral(@(t) exp(1i*n*t).*(K11(zk,nu,R,R2,t)-1i*sign(n)*K11H(zk,nu,R,R2,t)),0,2*pi)];
        k12s = [k12s; integral(@(t) exp(1i*n*t).*K12(zk,nu,R,R2,t),0,2*pi)];
        k21s = [k21s; integral(@(t) exp(1i*n*t).*(K21(zk,nu,R,R2,t)-1i*sign(n)*K21H(zk,nu,R,R2,t)),0,2*pi, "AbsTol",1E-8,"RelTol",1E-8)];
        k22s = [k22s; integral(@(t) exp(1i*n*t).*K22(zk,nu,R,R2,t),0,2*pi)];
    end
    
    k11 = pols.'*v2c*k11s;
    k12 = pols.'*v2c*k12s;
    k21 = pols.'*v2c*k21s;
    k22 = pols.'*v2c*k22s;
    
    lhs = [k11, k12; k21, k22];
    
    bc1 = integral(@(t) exp(1i*n*t).*BC1(zk,nu,sqrt(src(1)^2+src(2)^2),R,t).',0,2*pi);
    bc2 = integral(@(t) exp(1i*n*t).*BC2(zk,nu,sqrt(src(1)^2+src(2)^2),R,t).',0,2*pi);
    rhs = exp(-1i*n*atan2(src(2),src(1)))*[bc1;bc2];
    
    sol = lhs \ rhs;
    
    rho1 = rho1+sol(1)*exp(1i*n*thetas);
    Hrho1 = Hrho1-1i*sign(n)*sol(1)*exp(1i*n*thetas);
    rho2 = rho2+sol(2)*exp(1i*n*thetas);

end

srcs = R*[cos(thetas); sin(thetas)];

[val,grad] = flex2d.hkdiffgreen(zk,srcs,targs);

[~,nt] = size(targs);

nx = cos(thetas).';
ny = sin(thetas).';

taux = -sin(thetas).';
tauy = cos(thetas).';

nx = repmat(nx.',[nt 1]);
ny = repmat(ny.',[nt 1]);
taux = repmat(taux.',[nt 1]);
tauy = repmat(tauy.',[nt 1]);

K1 = -1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny);  % G_{ny}
K1H = -1/(2*zk^2).*(grad(:, :, 1).*taux + grad(:, :, 2).*tauy);  % G_{tauy}
K2 = 1/(2*zk^2).*val; % G

u(out) = K1*rho1.'*dt/(2*pi) + beta*K1H*Hrho1.'*dt/(2*pi) + K2*rho2.'*dt/(2*pi);

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

function val = K11(zk,nu,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,~,third] = flex2d.hkdiffgreen(zk,src,targ);

    nx = cos(t);
    ny = sin(t);

    nxtarg = 1;
    nytarg = 0;

    tauxtarg = 0;
    tauytarg = 1;
    
   val = -(1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny))) - ...
       nu./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ;  % first kernel with no hilbert transforms (G_{nx nx ny + nu G_{taux taux ny}).

end

function val = K11H(zk,nu,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,~,third] = flex2d.hkdiffgreen(zk,src,targ);

    nxtarg = 1;
    nytarg = 0;

    tauxtarg = 0;
    tauytarg = 1;

    taux = -sin(t);
    tauy = cos(t);
    
    val = -(1+ nu)/2*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) )  - ...
       (1+ nu)/2*nu.*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy))) ; % first kernel part coupled to hilbert transform (G_{nx nx tauy + nu G_{taux taux tauy}).
end


function val = K12(zk,nu,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,hess] = flex2d.hkdiffgreen(zk,src,targ);

    nx = 1;
    ny = 0;

    taux = 0;
    tauy = 1;
    
    val = 1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx+nu*taux.*taux) + ...
        hess(:, :, 2).*2.*(nx.*ny + nu*taux.*tauy)+ ...
        hess(:, :, 3).*(ny.*ny+nu*tauy.*tauy) );
end

function val = K21(zk,nu,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,~,third,forth] = flex2d.hkdiffgreen(zk,src,targ);

    nx = cos(t);
    ny = sin(t);

    nxtarg = 1;
    nytarg = 0;

    tauxtarg = 0;
    tauytarg = 1;
    
   val = 1/rs./(2*zk^2).*(1-nu).*((third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
            third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
           (third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ) ...
        - 1/(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          forth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) - ...
          ((2-nu)/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) ) ;

end

function val = K21H(zk,nu,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,~,third,forth] = flex2d.hkdiffgreen(zk,src,targ);

    taux = -sin(t);
    tauy = cos(t);

    nxtarg = 1;
    nytarg = 0;

    tauxtarg = 0;
    tauytarg = 1;
    
   val = 1/rs.*(1-nu).*(-((1+ nu)/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ nu)/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)))) ...
        -(1+ nu)/2.*(1/(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*taux) + ...
          forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*tauy + 3*nxtarg.*nxtarg.*nytarg.*taux) + ...
          forth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*tauy + 3*nxtarg.*nytarg.*nytarg.*taux) +...
          forth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*tauy + nytarg.*nytarg.*nytarg.*taux) +...
          forth(:, :, 5).*(nytarg.*nytarg.*nytarg.*tauy)) ) - ...
          ((2-nu)/2)*(1+nu).*(1/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*taux) + ...
          forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*tauy + tauxtarg.*tauxtarg.*nytarg.*taux + 2*tauxtarg.*tauytarg.*nxtarg.*taux) + ...
          forth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*tauy + 2*tauxtarg.*tauytarg.*nxtarg.*tauy + tauytarg.*tauytarg.*nxtarg.*taux + 2*tauxtarg.*tauytarg.*nytarg.*taux) +...
          forth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*tauy + 2*tauxtarg.*tauytarg.*nytarg.*tauy + tauytarg.*tauytarg.*nytarg.*taux) +...
         forth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*tauy)) ) ;
end

function val = K22(zk,nu,R1,rs,t)
    
    src = R1*[cos(t); sin(t)];
    targ = [rs 0*rs].';
    [~,~,hess,third] = flex2d.hkdiffgreen(zk,src,targ);

    nx = 1;
    ny = 0;

    taux = 0;
    tauy = 1;
    
    val = 1/(2*zk^2).*(third(:,:,1).*nx.^3 + third(:,:,2).*3.*nx.^2.*ny + ...
        + third(:,:,3).*3.*nx.*ny.^2 + third(:,:,4).*ny.^3) + ... 
        (2-nu)/(2*zk^2).*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(2.*nx.*taux.*tauy + ny.*taux.^2) + ...
        + third(:,:,3).*(2*ny.*taux.*tauy + nx.*tauy.^2) + third(:,:,4).*ny.*tauy.^2) + ...
        (1-nu)/rs/(2*zk^2).*(hess(:, :, 1).*(taux.*taux - nx.*nx) + ...
        hess(:, :, 2).*2.*(taux.*tauy-nx.*ny)+ ...
        hess(:, :, 3).*(tauy.*tauy-ny.*ny) );

end


function val = BC1(zk,nu,R1,rs,t)
    
    src = [R1 0*rs].';
    targ = rs*[cos(t); sin(t)];
    [~,~,hess] = flex2d.hkdiffgreen(zk,src,targ);

    nx = cos(t).';
    ny = sin(t).';

    taux = -sin(t).';
    tauy = cos(t).';
    
    val = 1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx+nu*taux.*taux) + ...
        hess(:, :, 2).*2.*(nx.*ny + nu*taux.*tauy)+ ...
        hess(:, :, 3).*(ny.*ny+nu*tauy.*tauy) );


% firstbc = 1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))-...
%            1/(2*zk^2).*(hessK(:, :, 1).*(nx.*nx) + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*(ny.*ny))+...
%            coefs(1)/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))-...
%            coefs(1)/(2*zk^2).*(hessK(:, :, 1).*(taux.*taux) + hessK(:, :, 2).*(2*taux.*tauy) + ...
%            hessK(:, :, 3).*(tauy.*tauy));

end

function val = BC2(zk,nu,R1,rs,t)
    
    src = [R1 0*rs].';
    targ = rs*[cos(t); sin(t)];
    [~,~,hess,third] = flex2d.hkdiffgreen(zk,src,targ);

    nx = cos(t).';
    ny = sin(t).';

    taux = -sin(t).';
    tauy = cos(t).';
    
    val = 1/(2*zk^2).*(third(:,:,1).*nx.^3 + third(:,:,2).*3.*nx.^2.*ny + ...
        + third(:,:,3).*3.*nx.*ny.^2 + third(:,:,4).*ny.^3) + ... 
        (2-nu)/(2*zk^2).*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(2.*nx.*taux.*tauy + ny.*taux.^2) + ...
        + third(:,:,3).*(2*ny.*taux.*tauy + nx.*tauy.^2) + third(:,:,4).*ny.*tauy.^2) + ...
        (1-nu)/rs/(2*zk^2).*(hess(:, :, 1).*(taux.*taux - nx.*nx) + ...
        hess(:, :, 2).*2.*(taux.*tauy-nx.*ny)+ ...
        hess(:, :, 3).*(tauy.*tauy-ny.*ny) );

% secondbc = 1./(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
% third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) - ...
%  1./(2*zk^2).*(thirdK(:, :, 1).*(nx.*nx.*nx) + thirdK(:, :, 2).*(3*nx.*nx.*ny)+...
%  thirdK(:, :, 3).*(3*nx.*ny.*ny) + thirdK(:, :, 4).*(ny.*ny.*ny)) +...
%  (2-coefs(1))/(2*zk^2).*(third(:, :, 1).*(taux.*taux.*nx) + third(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
%  third(:, :, 3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
%  + third(:, :, 4).*(tauy.*tauy.*ny)) - ...
%  (2-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(taux.*taux.*nx) + thirdK(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
%  thirdK(:, :, 3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
%  + thirdK(:, :, 4).*(tauy.*tauy.*ny)) +...
%  (1-coefs(1)).*(kappa).*(1/(2*zk^2).*(hess(:, :, 1).*taux.*taux + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*tauy.*tauy)-...
%   1/(2*zk^2).*(hessK(:, :, 1).*taux.*taux + hessK(:, :, 2).*(2*taux.*tauy) + ...
%  hessK(:, :, 3).*tauy.*tauy)-...
%  (1/(2*zk^2).*(hess(:, :, 1).*nx.*nx + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*ny.*ny)-...
% 1/(2*zk^2).*(hessK(:, :, 1).*nx.*nx + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*ny.*ny)));


end