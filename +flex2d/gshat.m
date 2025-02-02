function [val] = gshat(n,rho,rhop,rts,ejs,ifr2logr)
%
% computes the flexural gravity green's function for the 
% integro-differential equation determined by the roots of the polynomial:
%             z^5 - beta*z + gamma = 0
%
% output is a cell array with all the following kernels:
% - grad(:,:,1) has G_{x1}, grad(:,:,2) has G_{x2}
% - hess(:,:,1) has G_{x1x1}, hess(:,:,2) has G_{x1x2}, 
% hess(:,:,3) has G_{x2x2}
% - der3 has the third derivatives in the order G_{x1x1x1}, G_{x1x1x2}, 
% G_{x1x2x2}, G_{x2x2x2}
% - der4 has the fourth derivatives in the order G_{x1x1x1x1}, 
% G_{x1x1x1x2}, G_{x1x1x2x2}, G_{x1x2x2x2}, G_{x2x2x2x2}
%
% input:
%
% src - (2,ns) array of source locations
% targ - (2,nt) array of target locations
% rts - (5,) roots of polynomials above 
% ejs - (5,) coefficients of partial fraction expansion
%
% optional input:
%
% ifr2logr - boolean, default: false. If true, also subtract off the 
%             k^2/(8pi) r^2 log r kernel

if nargin < 6
    ifr2logr = false;
end

[~,ns] = size(rhop);
[~,nt] = size(rho);

rhop = repmat(rhop,nt,1);
rho = repmat(rho.',1,ns);
deltarho = rhop - rho;

np = 100;
thetas = -pi:(2*pi/np):(pi - pi/np);

val = 0*rho(:);
src = [0;0];

targs = sqrt(deltarho(:).^2 + 4*rho(:).*(rho(:) + deltarho(:))*sin(thetas/2).^2);

targs = [sqrt(deltarho(:).^2 + 4*rho(:).*(rho(:) + deltarho(:))*sin(thetas(k))^2)  rho(:)*0].';

val = val + exp(1i*n*thetas(k)).*flex2d.gs(src,targs,rts,ejs,ifr2logr);

val = pi*val.*rhop(:)/(2*np);

val = reshape(val,[nt ns]);

end