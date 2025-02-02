function [val,grad] = green(k,src,targ)
%FLEX2D.helm1d evaluate the helmholtz green's function in 1d
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

val = exp(1i*k*abs(rx))/(2i*k);

grad = 1/2*sign(rx).*exp(1i*k*abs(rx));
