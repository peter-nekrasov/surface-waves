function [xlegs, wlegs] = get_grid(zk, rt, dr, dz, ppw)

    if (nargin <= 4); ppw = 10; end
    persistent xlegloc wlegloc
    k = 16;
    if isempty(xlegloc) && isempty(wlegloc)
        [xlegloc,wlegloc] = lege.exps(k);
    end

    rs = rt + dr;
    drdiff = sqrt((rt+rs)^2  + dz^2) - sqrt(dr^2 + dz^2);
    npan = max(ceil(abs(zk)*drdiff*ppw/2/pi/k),3);
    
    h = pi/npan;
    tends = h:h:pi;
    
    dr0 = sqrt(dr^2 + dz^2);
    nref = max(ceil(-log(dr0/h)/log(2)),2);
    tends = [2.^(-nref:-1)*h tends];
    tstarts = [0 tends(1:end-1)];
    
    xlegs = tstarts + (xlegloc+1)/2*(tends-tstarts);
    wlegs = wlegloc/2*(tends-tstarts);
    xlegs = xlegs(:);
    wlegs = wlegs(:);
end
