function [val,grad] = gsaxisym(rts,ejs,n,src,targ)
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

[~,nt] = size(targ.r);
rhop = repmat(src.r(1,:),nt,1);

val = integral(@(t) exp(1i*n*t).*get_val(rts,ejs,src,targ,t),0,2*pi,'ArrayValued',true);
val = val.*rhop;

if nargout > 1

    grad = integral(@(t) exp(1i*n*t).*get_grad(rts,ejs,src,targ,t),0,2*pi,'ArrayValued',true);
    grad = grad.*rhop;

end

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


