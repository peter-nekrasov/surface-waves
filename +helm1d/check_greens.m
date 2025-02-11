% finite difference

h = 0.01;

src = [];
src.r = [5;0];

targ = [];
targ.r = [0:h:8*h; (0:h:8*h)*0];

d2 = [-1/560 8/315 -1/5 8/5 -205/72 8/5 -1/5 8/315 -1/560]/h^2;

% test on helmholtz green function first

k = 2;
[g,~] = helm1d.green(k,src,targ);

err = d2*g + k^2*g(5); % passed

% now check nonlocal greens function

beta = 2+0.4i;
gamma = 1;
[rts,ejs] = helm1d.find_roots(beta,gamma);

[gs,~] = helm1d.gshelm(rts,ejs,src,targ);
[gphi,~] = helm1d.gphihelm(rts,ejs,src,targ);

err = d2*gs + beta*gs(5) + gamma*gphi(5); % passed