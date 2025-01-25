function out = green(src,targ,rts,ejs)
%
% computes the green's function centered at (x,y) = 0 for the 
% integro-differential equation determined by the roots of the polynomial:
%             z^5 - beta*z + gamma = 0
%
% output is a cell array with all the kernels needed to solve the adjointed
% Lippman Schwinger equation: {val,hessxx,hessxy,hessyy,gradlapx,gradlapy}
% where:
% - val is the value of the Green's function centered at zero and
%   evaluated at (x,y)
% - hessxx is G_{xx}, hessxy is G_{xy}, 
%   hessyy is G_{yy}
% - gradlap is the gradient of the Laplacian, namely 
%   gradlapx is G_{xxx} + G_{xyy}, gradlapy is G_{yxx} + G_{yyy}
%
% input:
%
% x - x-coordinates array
% y - y-coordinates array
% beta - coefficient beta in the equation
% gamma - coefficient gamma in the equation
%
% optional input:
%
% opt - bool, default: false. 
%         Possible options are:
%         opt = false => Green's function (and derivatives)
%         opt = true => kernel used to evaluate phi on surface
%

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
% gradx = 0;
% grady = 0;
hessxx = 0;
hessxy = 0;
hessyy = 0;
gradlapx = 0;
gradlapy = 0;

for i = 1:5
    
    rhoj = rts(i);
    ej = ejs(i);

    if (angle(rhoj) == 0) && (rhoj ~= 0)

       [sk0,~,hesssk0,gradlapsk0] = struveKdiffgreen(rhoj,src,targ);
       [h0,~,hessh0,thirdh0] = helmdiffgreen(rhoj,src,targ);

       h0(r == 0) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));

       h0 = -4i*h0;
       % gradh0 = -4i*gradh0;
       
       % h0x = gradh0(:,:,1);
       % h0y = gradh0(:,:,2);
       % 
       % h0x(r == 0) = 0;
       % h0y(r == 0) = 0;
       
       h0xx = hessh0(:,:,1);
       h0xy = hessh0(:,:,2);
       h0yy = hessh0(:,:,3);

       h0xx(r == 0) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
       h0xy(r == 0) = 0;
       h0yy(r == 0) = rhoj^2/(4*pi)*(log(rhoj)-1i*pi/2+eulergamma-0.5-log(2));
       
       h0xx = -4i*h0xx;
       h0xy = -4i*h0xy;
       h0yy = -4i*h0yy;

       h0xxx = thirdh0(:,:,1);
       h0yxx = thirdh0(:,:,2);
       h0xyy = thirdh0(:,:,3);
       h0yyy = thirdh0(:,:,4);

       h0xxx(r == 0) = 0;
       h0yxx(r == 0) = 0;
       h0xyy(r == 0) = 0;
       h0yyy(r == 0) = 0;

       h0xxx = -4i*h0xxx;
       h0yxx = -4i*h0yxx;
       h0xyy = -4i*h0xyy;
       h0yyy = -4i*h0yyy;
       
       % sk0x = gradsk0(:,:,1);
       % sk0y = gradsk0(:,:,2);

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);

       sk0lapx = gradlapsk0(:,:,1);
       sk0lapy = gradlapsk0(:,:,2);

       val = val + ej*rhoj^2*(-sk0 + 2i*h0);
       phi = phi + ej*rhoj*(-sk0 + 2i*h0);

       %gradx = gradx + ej*rhoj^2*(-sk0x + 2i*h0x);
       %grady = grady + ej*rhoj^2*(-sk0y + 2i*h0y);
       
       hessxx = hessxx + ej*rhoj^2*(-sk0xx + 2i*h0xx);
       hessxy = hessxy + ej*rhoj^2*(-sk0xy + 2i*h0xy);
       hessyy = hessyy + ej*rhoj^2*(-sk0yy + 2i*h0yy);

       gradlapx = gradlapx + ej*rhoj^2*(-sk0lapx + 2i*(h0xxx+h0xyy));
       gradlapy = gradlapy + ej*rhoj^2*(-sk0lapy + 2i*(h0yxx+h0yyy));

    elseif rhoj ~= 0

       [sk0,~,hesssk0,gradlapsk0] = struveKdiffgreen(-rhoj,src,targ);

       % sk0x = gradsk0(:,:,1);
       % sk0y = gradsk0(:,:,2);

       sk0xx = hesssk0(:,:,1);
       sk0xy = hesssk0(:,:,2);
       sk0yy = hesssk0(:,:,3);

       sk0lapx = gradlapsk0(:,:,1);
       sk0lapy = gradlapsk0(:,:,2);

       val = val + ej*rhoj^2*sk0;
       phi = phi + ej*rhoj*sk0;

       %gradx = gradx + ej*rhoj^2*sk0x;
       %grady = grady + ej*rhoj^2*sk0y;

       hessxx = hessxx + ej*rhoj^2*sk0xx;
       hessxy = hessxy + ej*rhoj^2*sk0xy;
       hessyy = hessyy + ej*rhoj^2*sk0yy;

       gradlapx = gradlapx + ej*rhoj^2*sk0lapx;
       gradlapy = gradlapy + ej*rhoj^2*sk0lapy;

    end

end

val = 1/2*val;
phi = 1/2*phi;
%gradx = 1/2*gradx;
%grady = pi/2*grady;
hessxx = 1/2*hessxx;
hessxy = 1/2*hessxy;
hessyy = 1/2*hessyy;
gradlapx = 1/2*gradlapx;
gradlapy = 1/2*gradlapy;

hess = cat(3,hessxx,hessxy,hessyy);
gradlap = cat(3,gradlapx,gradlapy);

out = {val,hess,gradlap,phi};

end