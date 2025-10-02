function [val,grad,lap] = gsaxisym(rts,ejs,n,src,targ)
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
% src - source array - rho' in first row, zeros in second row
% targ - target array - rho in first row, zeros in second row
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

[val,grad,hess] = helm2d.gshelm(rts,ejs,src,tar);

eintw = cos(n*xlegs).*wlegs;

val = reshape(val.',[nt*ns ntheta]);
val = 2*val*eintw;
val = reshape(val,[ns nt]).';
val = val.*rhop;

if nargout > 1

    grad = grad(:,:,1).*cos(thetas(:))+grad(:,:,2).*sin(thetas(:));
    grad = reshape(grad.',[nt*ns ntheta]);
    grad = 2*grad*eintw;
    grad = reshape(grad,[ns nt]).';
    grad = grad.*rhop;

end
if nargout > 2

    lap = reshape(hess(:,:,1).'+hess(:,:,3).',[nt*ns ntheta]);
    lap = 2*lap*eintw;
    lap = reshape(lap,[ns nt]).';
    lap = lap.*rhop;

end

% val2 = integral(@(t) exp(1i*n*t).*get_val(rts,ejs,src,targ,t),0,2*pi,'ArrayValued',true);
% val2 = val2.*rhop;
% 
% grad2 = integral(@(t) exp(1i*n*t).*get_grad(rts,ejs,src,targ,t),0,2*pi,'ArrayValued',true);
% grad2 = grad2.*rhop;
% 
% lap2 = integral(@(t) exp(1i*n*t).*get_lap(rts,ejs,src,targ,t),0,2*pi,'ArrayValued',true);
% lap2 = lap2.*rhop;
% 
% nnn = 3;

end


function val = get_val(rts,ejs,src,targ,theta)

    targ.r = [targ.r(1,:)*cos(theta); targ.r(1,:)*sin(theta)];
    val = helm2d.gshelm(rts,ejs,src,targ);

end

function grad = get_grad(rts,ejs,src,targ,theta)

    targ.r = [targ.r(1,:)*cos(theta); targ.r(1,:)*sin(theta)];
    [~,grad] = helm2d.gshelm(rts,ejs,src,targ);
    grad = grad(:,:,1)*cos(theta) + grad(:,:,2)*sin(theta);

end

function lap = get_lap(rts,ejs,src,targ,theta)

    targ.r = [targ.r(1,:)*cos(theta); targ.r(1,:)*sin(theta)];
    [~,~,hess] = helm2d.gshelm(rts,ejs,src,targ);
    lap = hess(:,:,1) + hess(:,:,3);

end