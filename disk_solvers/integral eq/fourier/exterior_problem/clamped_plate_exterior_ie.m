clear
close all 
addpath(genpath('../..'))

[rs,~,v2c] = lege.exps(16);
[pols,~] = lege.pols(-1,15);

rs = rs/8 + 1/8;

zk = 2;
R = sqrt(2);


% calculate fourier modes 

dt = 0.001;
thetas = dt:dt:2*pi;
rho1 = thetas*0;
rho2 = thetas*0;

N = 10;

src = [0;0];
targ = [R;0];

% [val,grad] = flex2d.hkdiffgreen(zk,src,targ);
% 
% bc1 = 1/(2*zk^2)*val;
% bc2 = 1/(2*zk^2)*grad(1,1,1);

for n = -N:N

    k11s = [];
    k12s = [];
    k21s = [];
    k22s = [];
    
    for i = 1:16
        R2 = R+rs(i);
        k11s = [k11s; integral(@(t) R*exp(1i*n*t).*K11(zk,R,R2,t),0,2*pi)];
        k12s = [k12s; integral(@(t) R*exp(1i*n*t).*K12(zk,R,R2,t),0,2*pi)];
        k21s = [k21s; integral(@(t) R*exp(1i*n*t).*K21(zk,R,R2,t),0,2*pi)];
        k22s = [k22s; integral(@(t) R*exp(1i*n*t).*K22(zk,R,R2,t),0,2*pi)];
    end
    
    k11 = pols.'*v2c*k11s;
    k12 = pols.'*v2c*k12s;
    k21 = pols.'*v2c*k21s;
    k22 = pols.'*v2c*k22s;
    
    lhs = [k11, k12; k21, k22];

    [val,grad] = flex2d.hkdiffgreen(zk,src,targ);
    
    bc1 = integral(@(t) 1/(2*pi)*exp(1i*n*t).*BC1(zk,sqrt(src(1)^2+src(2)^2),R,t),0,2*pi);
    bc2 = integral(@(t) 1/(2*pi)*exp(1i*n*t).*BC2(zk,sqrt(src(1)^2+src(2)^2),R,t),0,2*pi);
    rhs = exp(-1i*n*atan2(src(2),src(1)))*[bc1;bc2];
    
    sol = lhs \ rhs;
    

    rho1 = rho1+sol(1)*exp(1i*n*thetas);
    rho2 = rho2+sol(2)*exp(1i*n*thetas);

end

new_targs = [R:0.1:R+1;R:0.1:R+1];
srcs = R*[cos(thetas); sin(thetas)];

[~,~,hess,third] = flex2d.hkdiffgreen(zk,srcs,new_targs);

[~,nt] = size(new_targs);

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

utarg = K1*rho1.'*R*dt + K2*rho2.'*R*dt;

[true_sol,~] = flex2d.hkdiffgreen(zk,src,new_targs);
true_sol = true_sol/(2*zk^2);

err = max(abs(utarg - true_sol)) ./ max(abs(true_sol))


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