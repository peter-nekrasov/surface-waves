function [val] = modal_gf(k,src,targ,n)
%MODAL_GF evaluate the modal Green's function, i.e. 
%
% G(x,y) = pi*rho'/(2*k^2)*(i/4 H_n^(1)(k*rho>) J_n(k*rho<) -
%                                        1/(2 pi) K_n(k rho>) I_n(k rho<)))
%
% where rho> is the larger radius, rho< is the lesser radius and
% src = [rho,z], targ = [rho,z]
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
%
% optional input:
%
% ifr2logr - boolean, default: false. If true, also subtract off the 
%             k^2/(8pi) r^2 log r kernel

if nargin < 4
    n = 0;
end

src = src.r;
targ = targ.r;

[~,ns] = size(src);
[~,nt] = size(targ);

rhos = repmat(src(1,:),nt,1);
rhot = repmat(targ(1,:).',1,ns);

rhogreater = rhos;
rholesser = rhos;

rhogreater(rhot > rhos) = rhot(rhot > rhos);
rholesser(rhot <= rhos) = rhot(rhot <= rhos);

val = rhos.*pi./k^2.*(1i/4*besselh(n,k*rhogreater).*besselj(n,k*rholesser) ...
    - 1/(2*pi)*besselk(n,k*rhogreater).*besseli(n,k*rholesser));

end

