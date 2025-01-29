function [val,grad,hess,third,fourth] = fggreen(src,targ,rts,ejs,ifr2logr)
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

if nargin < 5
    ifr2logr = false;
end


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

val = 0;
phi = 0;

gradx = 0;
grady = 0;

hessxx = 0;
hessxy = 0;
hessyy = 0;

thirdxxx = 0;
thirdxxy = 0;
thirdxyy = 0;
thirdyyy = 0;

fourthxxxx = 0;
fourthxxxy = 0;
fourthxxyy = 0;
fourthxyyy = 0;
fourthyyyy = 0;


for i = 1:5
    
    rhoj = rts(i);
    ej = ejs(i);

    if (angle(rhoj) == 0) && (rhoj ~= 0)

       [sk0,gradsk0,hesssk0,thirdsk0,fourthsk0] = flex2d.struveKdiffgreen(rhoj,src,targ,ifr2logr);
       [h0,gradh0,hessh0,thirdh0,fourthh0] = flex2d.helmdiffgreen(rhoj,src,targ,ifr2logr);

       h0(r == 0) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));

       h0 = -4i*h0;
       gradh0 = -4i*gradh0;
       
       h0x = gradh0(:,:,1);
       h0y = gradh0(:,:,2);
       
       h0xx = hessh0(:,:,1);
       h0xy = hessh0(:,:,2);
       h0yy = hessh0(:,:,3);

       h0xxx = thirdh0(:,:,1);
       h0xxy = thirdh0(:,:,2);
       h0xyy = thirdh0(:,:,3);
       h0yyy = thirdh0(:,:,4);

       h0xxxx = fourthh0(:,:,1);
       h0xxxy = fourthh0(:,:,2);
       h0xxyy = fourthh0(:,:,3);
       h0xyyy = fourthh0(:,:,4);
       h0yyyy = fourthh0(:,:,5);

       h0x(r == 0) = 0;
       h0y(r == 0) = 0;

       h0xx(r == 0) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
       h0xy(r == 0) = 0;
       h0yy(r == 0) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
       
       h0xx = -4i*h0xx;
       h0xy = -4i*h0xy;
       h0yy = -4i*h0yy;

       h0xxx(r == 0) = 0;
       h0xxy(r == 0) = 0;
       h0xyy(r == 0) = 0;
       h0yyy(r == 0) = 0;

       h0xxx = -4i*h0xxx;
       h0xxy = -4i*h0xxy;
       h0xyy = -4i*h0xyy;
       h0yyy = -4i*h0yyy;
       
       h0xxxx = -4i*h0xxxx;
       h0xxxy = -4i*h0xxxy;
       h0xxyy = -4i*h0xxyy;
       h0xyyy = -4i*h0xyyy;
       h0yyyy = -4i*h0yyyy;

       sk0x = gradsk0(:,:,1);
       sk0y = gradsk0(:,:,2);

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);

       sk0xxx = thirdsk0(:,:,1);
       sk0xxy = thirdsk0(:,:,2);
       sk0xyy = thirdsk0(:,:,3);
       sk0yyy = thirdsk0(:,:,4);

       sk0xxxx = fourthsk0(:,:,1);
       sk0xxxy = fourthsk0(:,:,2);
       sk0xxyy = fourthsk0(:,:,3);
       sk0xyyy = fourthsk0(:,:,4);
       sk0yyyy = fourthsk0(:,:,5);

       val = val + ej*rhoj^2*(-sk0 + 2i*h0);
       phi = phi + ej*rhoj*(-sk0 + 2i*h0);

       gradx = gradx + ej*rhoj^2*(-sk0x + 2i*h0x);
       grady = grady + ej*rhoj^2*(-sk0y + 2i*h0y);
       
       hessxx = hessxx + ej*rhoj^2*(-sk0xx + 2i*h0xx);
       hessxy = hessxy + ej*rhoj^2*(-sk0xy + 2i*h0xy);
       hessyy = hessyy + ej*rhoj^2*(-sk0yy + 2i*h0yy);

       thirdxxx = thirdxxx + ej*rhoj^2*(-sk0xxx + 2i*h0xxx);
       thirdxxy = thirdxxy + ej*rhoj^2*(-sk0xxy + 2i*h0xxy);
       thirdxyy = thirdxyy + ej*rhoj^2*(-sk0xyy + 2i*h0xyy);
       thirdyyy = thirdyyy + ej*rhoj^2*(-sk0yyy + 2i*h0yyy);

       fourthxxxx = fourthxxxx + ej*rhoj^2*(-sk0xxxx + 2i*h0xxxx);
       fourthxxxy = fourthxxxy + ej*rhoj^2*(-sk0xxxy + 2i*h0xxxy);
       fourthxxyy = fourthxxyy + ej*rhoj^2*(-sk0xxyy + 2i*h0xxyy);
       fourthxyyy = fourthxyyy + ej*rhoj^2*(-sk0xyyy + 2i*h0xyyy);
       fourthyyyy = fourthyyyy + ej*rhoj^2*(-sk0yyyy + 2i*h0yyyy);

    elseif rhoj ~= 0

       [sk0,gradsk0,hesssk0,thirdsk0,fourthsk0] = flex2d.struveKdiffgreen(-rhoj,src,targ,ifr2logr);

       sk0x = gradsk0(:,:,1);
       sk0y = gradsk0(:,:,2);

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);

       sk0xxx = thirdsk0(:,:,1);
       sk0xxy = thirdsk0(:,:,2);
       sk0xyy = thirdsk0(:,:,3);
       sk0yyy = thirdsk0(:,:,4);

       sk0xxxx = fourthsk0(:,:,1);
       sk0xxxy = fourthsk0(:,:,2);
       sk0xxyy = fourthsk0(:,:,3);
       sk0xyyy = fourthsk0(:,:,4);
       sk0yyyy = fourthsk0(:,:,5);

       val = val + ej*rhoj^2*sk0;
       phi = phi + ej*rhoj*sk0;

       gradx = gradx + ej*rhoj^2*sk0x;
       grady = grady + ej*rhoj^2*sk0y;

       hessxx = hessxx + ej*rhoj^2*sk0xx;
       hessxy = hessxy + ej*rhoj^2*sk0xy;
       hessyy = hessyy + ej*rhoj^2*sk0yy;

       thirdxxx = thirdxxx + ej*rhoj^2*sk0xxx;
       thirdxxy = thirdxxy + ej*rhoj^2*sk0xxy;
       thirdxyy = thirdxyy + ej*rhoj^2*sk0xyy;
       thirdyyy = thirdyyy + ej*rhoj^2*sk0yyy;

       fourthxxxx = fourthxxxx + ej*rhoj^2*sk0xxxx;
       fourthxxxy = fourthxxxy + ej*rhoj^2*sk0xxxy;
       fourthxxyy = fourthxxyy + ej*rhoj^2*sk0xxyy;
       fourthxyyy = fourthxyyy + ej*rhoj^2*sk0xyyy;
       fourthyyyy = fourthyyyy + ej*rhoj^2*sk0yyyy;

    end

end

val = 1/2*val;

grad(:,:,1) = 1/2*gradx;
grad(:,:,2) = 1/2*grady;

hess(:,:,1) = 1/2*hessxx;
hess(:,:,2) = 1/2*hessxy;
hess(:,:,3) = 1/2*hessyy;

third(:,:,1) = 1/2*thirdxxx;
third(:,:,2) = 1/2*thirdxxy;
third(:,:,3) = 1/2*thirdxyy;
third(:,:,4) = 1/2*thirdyyy;

fourth(:,:,1) = 1/2*fourthxxxx;
fourth(:,:,2) = 1/2*fourthxxxy;
fourth(:,:,3) = 1/2*fourthxxyy;
fourth(:,:,4) = 1/2*fourthxyyy;
fourth(:,:,5) = 1/2*fourthyyyy;

end