function [val] = ext_modal_gf(k,R,nu,src,targ,n)
%BDD_MODAL_GF evaluate the exterior modal Green's function, i.e. 
%
% G(rho,rho') = pi*rho'/(2*k^2)*(i/4 H_n^(1)(k*rho>) J_n(k*rho<) +
%          - 1/(2 pi) K_n(k rho>) I_n(k rho<))) 
%           + alpha(rho') H_n(k rho) + beta(rho') K_n(k rho) 
%
% where rho> is the larger radius, rho< is the lesser radius,
% src = [rho,z], targ = [rho,z], and alpha, beta are chosen to satisfy the 
% free plate boundary conditions at rho = R;
%
% where H_0^(1) is the principal branch of the Hankel function
% of the first kind and K_0 is the modified Bessel function of the second
% kind. This routine avoids numerical cancellation
% when |k||x-y| is small.
%
% - grad(:,:,1) has G_{x1}, grad(:,:,2) has G_{x2}
% - hess(:,:,1) has G_{x1x1}, hess(:,:,2) has G_{x1x2}, 
% hess(:,:,3) has G_{x2x2}
% - der3 has the third derivatives in the order G_{x1x1x1}, G_{x1x1x2}, 
% G_{x1x2x2}, G_{x2x2x2}
% - der4 has the fourth derivatives in the order G_{x1x1x1x1}, 
% G_{x1x1x1x2}, G_{x1x1x2x2}, G_{x1x2x2x2}, G_{x2x2x2x2}
%
% derivatives are on the *target* variables
%
% input:
%
% src - (2,ns) array of source locations
% targ - (2,nt) array of target locations
% k - wave number, as above
% n - mode number

if nargin < 6
    n = 0;
end

src = src.r;
targ = targ.r;

[~,ns] = size(src);
[~,nt] = size(targ);

rhop = repmat(src(1,:),nt,1);
rho = repmat(targ(1,:).',1,ns);

rhogreater = rhop;
rholesser = rhop;

rhogreater(rho > rhop) = rho(rho > rhop);
rholesser(rho <= rhop) = rho(rho <= rhop);

val = rhop.*pi./k^2.*(1i/4*besselh(n,k*rhogreater).*besselj(n,k*rholesser) ...
    - 1/(2*pi)*besselk(n,k*rhogreater).*besseli(n,k*rholesser));

hnR = besselh(n,k*R);
hnRp = k/2*(besselh(n-1,k*R)-besselh(n+1,k*R));
hnRpp = (-k*R*besselh(n-1,k*R)+(n+n^2-k^2*R^2)*besselh(n,k*R))/R^2;
hnRppp = (k*R*(2+n^2-k^2*R^2)*besselh(-1+n,k*R)...
           - (1+n)*(2*n+n^2-k^2*R^2)*besselh(n,k*R))/R^3;

jnR = besselj(n,k*R);
jnRp = k/2*(besselj(n-1,k*R)-besselj(n+1,k*R));
jnRpp = (-k*R*besselj(n-1,k*R)+(n+n^2-k^2*R^2)*besselj(n,k*R))/R^2;
jnRppp = (k*R*(2+n^2-k^2*R^2)*besselj(-1+n,k*R)...
           - (1+n)*(2*n+n^2-k^2*R^2)*besselj(n,k*R))/R^3;

knR = besselk(n,k*R);
knRp = k/2*(-besselk(n-1,k*R)-besselk(n+1,k*R));
knRpp = (k*R*besselk(n-1,k*R)+(n+n^2+k^2*R^2)*besselk(n,k*R))/R^2;
knRppp = (-k*R*(2+n^2+k^2*R^2)*besselk(-1+n,k*R)...
            -(1+n)*(2*n+n^2+k^2*R^2)*besselk(n,k*R))/R^3;

inR = besseli(n,k*R);
inRp = k/2*(besseli(n-1,k*R)+besseli(n+1,k*R));
inRpp = -(k*R*besseli(n-1,k*R)-2*(n^2+k^2*R^2)*besseli(n,k*R)+ k*R*besseli(n+1,k*R))/(2*R^2);
inRppp = (k*R*(2+n^2+k^2*R^2)*besseli(-1+n,k*R)...
            - 2*(3*n^2+k^2*R^2)*besseli(n, k*R)+ ...
            k*R*(2+n^2+k^2*R^2)*besseli(1 + n, k*R))/(2*R^3);

L1hnR = hnRpp + nu/R*hnRp - nu/R^2*n^2*hnR;
L1knR = knRpp + nu/R*knRp - nu/R^2*n^2*knR;
L1jnR = jnRpp + nu/R*jnRp - nu/R^2*n^2*jnR;
L1inR = inRpp + nu/R*inRp - nu/R^2*n^2*inR;
L2hnR = hnRppp + 1/R*hnRpp - (1+n^2*(2-nu))/R^2*hnRp + (3-nu)/R^3*n^2*hnR;
L2knR = knRppp + 1/R*knRpp - (1+n^2*(2-nu))/R^2*knRp + (3-nu)/R^3*n^2*knR;
L2jnR = jnRppp + 1/R*jnRpp - (1+n^2*(2-nu))/R^2*jnRp + (3-nu)/R^3*n^2*jnR;
L2inR = inRppp + 1/R*inRpp - (1+n^2*(2-nu))/R^2*inRp + (3-nu)/R^3*n^2*inR;

hnrhop = besselh(n,k*rhop);
knrhop = besselk(n,k*rhop);

L1GnR = pi*rhop/k^2.*(1i/4*hnrhop*L1jnR - 1/(2*pi)*knrhop*L1inR);
L2GnR = pi*rhop/k^2.*(1i/4*hnrhop*L2jnR - 1/(2*pi)*knrhop*L2inR);

alpha = 1/(L1knR*L2hnR - L2knR*L1hnR)*(L2knR*L1GnR - L1knR*L2GnR);
beta = 1/(L1knR*L2hnR - L2knR*L1hnR)*(-L2hnR*L1GnR + L1hnR*L2GnR);

val = val + alpha.*besselh(n,k*rho) + beta.*besselk(n,k*rho);


end

