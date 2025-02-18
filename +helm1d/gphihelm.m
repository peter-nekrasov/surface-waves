function [val,grad,hess] = gphihelm(rts,ejs,src,targ)
%HELM1d.gshelm evaluate the nonlocal helmholtz green's function in 1d
% for the given sources and targets
%
% All derivatives are w.r.t. the target xt

src = src.r(1,:);
targ = targ.r(1,:);

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src,nt,1);
xt = repmat(targ.',1,ns);

rx = xt-xs;

val = 0;
grad = 0;
hess = 0;

for i = 1:3

    rhoj = rts(i);
    ej = ejs(i);

    val = val + 1/(2*pi)*ej*(helm1d.expeval(1i*rhoj*rx)+helm1d.expeval(-1i*rhoj*rx));
    grad = grad + 1i*sign(rx)/(2*pi)*ej*rhoj.*(helm1d.expeval(1i*rhoj*abs(rx))-helm1d.expeval(-1i*rhoj*abs(rx)));
    hess = hess - 1/(2*pi)*ej*rhoj^2.*(helm1d.expeval(1i*rhoj*abs(rx))+helm1d.expeval(-1i*rhoj*abs(rx)));

    if (angle(rhoj) < pi/2) && (angle(rhoj) >= 0) && (rhoj ~= 0)
        
        val = val + 1i*ej*exp(1i*rhoj*abs(rx));
        grad = grad - ej*rhoj*sign(rx).*exp(1i*rhoj*abs(rx));    
        hess = hess - 1i*ej*rhoj^2*exp(1i*rhoj*abs(rx));

    elseif (angle(rhoj) > -pi/2) && (angle(rhoj) < 0) && (rhoj ~= 0)

        val = val + 1i*ej*exp(-1i*rhoj*abs(rx));
        grad = grad + ej*rhoj*sign(rx).*exp(-1i*rhoj*abs(rx));
        hess = hess + 1i*ej*rhoj^2*exp(-1i*rhoj*abs(rx));

    end

    val(rx == 0) = val(rx == 0) - ej.*(log(1i*rhoj)+log(-1i*rhoj))/(2*pi); 

end

grad(rx == 0) = 0; %1/2*sign(rx).*exp(1i*k*abs(rx));

end