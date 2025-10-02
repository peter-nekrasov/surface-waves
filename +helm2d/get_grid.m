function [xlegs, wlegs] = get_grid()

    persistent xlegloc wlegloc
    k = 16;
    if isempty(xlegloc) && isempty(wlegloc)
        [xlegloc,wlegloc] = lege.exps(k);
    end

    npan = 4;
    
    h = pi/npan;
    tends = h:h:pi;
    
    nref = 15;
    tends = [2.^(-nref:-1)*h tends];
    tstarts = [0 tends(1:end-1)];
    
    xlegs = tstarts + (xlegloc+1)/2*(tends-tstarts);
    wlegs = wlegloc/2*(tends-tstarts);
    xlegs = xlegs(:);
    wlegs = wlegs(:);
end
