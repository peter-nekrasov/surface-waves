function [val,grad,hess,gradlap] = struveKdiffgreen(rhoj,src,targ)
% evaluator for struveK(0,rhoj*x) and its derivatives with singularity
% subtraction
%
% K(x,y) = -i*R_0(rhoj|x-y|) + i*H^{(1)}_0(rhoj|x-y|) + 2/pi log(|x-y|)
%
% where the Hankel part has log subtracted off
% 
% output : (note the convention is not the same as helmdiffgreen)
%
% - val has the value of the Green's function centered at zero and
%   evaluated at (x,y)
% - grad(:,:,1) has G_{x}, grad(:,:,2) has G_{y}
% - hess(:,:,1) has G_{xx}, hess(:,:,2) has G_{xy}, 
%   hess(:,:,3) has G_{yy}
% - gradlap is the gradient of the Laplacian, namely 
%   gradlap(:,:,1) has G_{xxx} + G_{xyy}, gradlap(:,:,2) has G_{yxx} + G_{yyy}
%
% input: 
%
% rhoj - complex number
% src - [2,n] array of source locations
% targ - [2,n] array of target locations
%
% Note: this code has only been tested for src = [0,0]

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

dx2 = dx.*dx;
dy2 = dy.*dy;

r2 = dx2 + dy2;
r = sqrt(r2);
r3 = r.^3;

ilow = (imag(rhoj) < 0);

if ilow
    rhoj = conj(rhoj); % flip rhoj up into upper half plane
end

zt = r*rhoj;
[cr0,cr1] = struveR(zt);

[h0,gradh0,hessh0,thirdh0] = helmdiffgreen(rhoj,src,targ);

h0x = gradh0(:,:,1);
h0y = gradh0(:,:,2);

h0xx = hessh0(:,:,1);
h0xy = hessh0(:,:,2);
h0yy = hessh0(:,:,3);

h0xxx = thirdh0(:,:,1);
h0yxx = thirdh0(:,:,2);
h0xyy = thirdh0(:,:,3);
h0yyy = thirdh0(:,:,4);

h0(zt == 0) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
h0 = -4i*h0;

h0x = -4i*h0x;
h0y = -4i*h0y;

h0xx = -4i*h0xx;
h0xy = -4i*h0xy;
h0yy = -4i*h0yy;

h0xxx = -4i*h0xxx;
h0yxx = -4i*h0yxx;
h0xyy = -4i*h0xyy;
h0yyy = -4i*h0yyy;

val = -1i*cr0+1i*h0;
gradx = 1i*rhoj*cr1.*dx./r+1i*h0x;
grady = 1i*rhoj*cr1.*dy./r+1i*h0y;
hessxx = 1i*rhoj^2*dx.^2./r2.*cr0-1i*rhoj*(dx.^2 - dy.^2)./r3.*cr1+1i*h0xx; 
hessxy = 1i*rhoj^2*dx.*dy./r2.*cr0-2i*rhoj*dx.*dy./r3.*cr1+1i*h0xy; 
hessyy = 1i*rhoj^2*dy.^2./r2.*cr0+1i*rhoj*(dx.^2 - dy.^2)./r3.*cr1+1i*h0yy; 
gradlapx = 1i*rhoj^3*dx./r.*(-cr1) + 1i*(h0xxx+h0xyy);
gradlapy = 1i*rhoj^3*dy./r.*(-cr1) + 1i*(h0yxx+h0yyy);

% tst = gradient(hessxx+hessyy,0.01);
% tiledlayout(2,1)
% nexttile
% plot(xt,real(tst),xt,imag(tst))
% 
% nexttile
% plot(xt,real(gradlapx),xt,imag(gradlapx))

gradx(r == 0) = 0;
grady(r == 0) = 0;

hessxx(r == 0) = 1i*rhoj^2-1i*rhoj^2/2-4i*1i*rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
hessxy(r == 0) = 0;
hessyy(r == 0) = 1i*rhoj^2-1i*rhoj^2/2-4i*1i*rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));

gradlapx(r == 0) = 0;
gradlapy(r == 0) = 0;

if ilow
    val = conj(val);
    gradx = conj(gradx);
    grady = conj(grady);
    hessxx = conj(hessxx);
    hessxy = conj(hessxy);
    hessyy = conj(hessyy);
    gradlapx = conj(gradlapx);
    gradlapy = conj(gradlapy);
end

grad = cat(3,gradx,grady);
hess = cat(3,hessxx,hessxy,hessyy);
gradlap = cat(3,gradlapx,gradlapy);

end