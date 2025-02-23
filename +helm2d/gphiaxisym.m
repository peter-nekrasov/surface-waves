function [val,grad,lap] = gphiaxisym(rts,ejs,n,src,targ)
%
% computes the green's function for the integro-differential equation 
% determined by the polynomial:
%             alpha z^2 + beta + gamma/|z| = -2
%
% outputs are:
% - val is the value of the Green's function centered at zero and
%   evaluated at (x,y)
% - grad(:,:,1) is G_{x}, grad(:,:,2) is G_{y} 
% - hess(:,:,1) is G_{xx}, hess(:,:,2) is G_{xy}, 
%   hess(:,:,3) is G_{yy}
%
% input:
%
% src - source array
% targ - target array
% rts - roots of cubic polynomial
% ejs - residues (see notes)


[~,ns] = size(src.r);
[~,nt] = size(targ.r);

rhop = repmat(src.r(1,:),nt,1);

[xlegs,wlegs] = helm2d.get_grid();

ntheta = numel(xlegs);

tarx = targ.r(1,:).'*cos(xlegs).';
tary = targ.r(1,:).'*sin(xlegs).';
thetas = ones(nt,1)*xlegs.';

tar.r = [tarx(:) tary(:)].';

[val,grad,hess] = helm2d.gphihelm(rts,ejs,src,tar);

eintw = cos(n*xlegs).*wlegs;

val = reshape(val.',[nt*ns ntheta]);
val = 2*val*eintw;
val = reshape(val,[nt ns]);
val = val.*rhop;

if nargout > 1

    grad = grad(:,:,1).*cos(thetas(:))+grad(:,:,2).*sin(thetas(:));
    grad = reshape(grad.',[nt*ns ntheta]);
    grad = 2*grad*eintw;
    grad = reshape(grad,[nt ns]);
    grad = grad.*rhop;

end
if nargout > 2

    lap = reshape(hess(:,:,1).'+hess(:,:,3).',[nt*ns ntheta]);
    lap = 2*lap*eintw;
    lap = reshape(lap,[nt ns]);
    lap = lap.*rhop;

end

end