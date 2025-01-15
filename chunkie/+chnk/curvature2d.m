function kappa = curvature2d(ptinfo)
%CHNK.CURVATURE2D signed curvature of points on 2D curve
%
%   kappa = (x'y''-x''y')/( sqrt(x'^2+y'^2)^3 )
% 
% Syntax: kappa = chnk.curvature2d(ptinfo)
%
% Input:
%   ptinfo - curve point info struct, with entries
%       ptinfo.r - positions (2,:) array
%       ptinfo.d - first derivative in underlying parameterization (2,:)
%       ptinfo.d2 - second derivative in underlying parameterization (2,:)
%
% Output:
%   nrm - (2,:) array containing corresponding normal information
%
% see also CHNK.NORMAL2D

d = ptinfo.d;
sh = size(d); sh = sh(2:end);
d2 = ptinfo.d2;
dnrm3 = sqrt(sum(d.^2,1)).^3;
kappa = bsxfun(@rdivide,d(1,:).*d2(2,:)-d(2,:).*d2(1,:),dnrm3(:).');
if (length(sh) > 1) kappa = reshape(kappa,sh); end