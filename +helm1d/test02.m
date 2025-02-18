%% test for derivatives of greens functions
% part 1 - check helmholtz kernels
% 

k = 3;

src = [];
src.r = [5;0];

h = 0.01;

targ = [];
targ.r = [0:h:8*h; (0:h:8*h)*0];

d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;
d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560	]/h^2;

[val,grad,hess] = helm1d.green(k,src,targ);
grad = grad(:,:,1);

err = d1*val - grad(5);
err2 = d2*val - hess(5);


%% part 2 - nonlocal helmholtz kernels (G_S)


beta = 1;
gamma = 2;
[rts,ejs] = helm1d.find_roots(beta,gamma);

src = [];
src.r = [5;0];

h = 0.001;

targ = [];
targ.r = [0:h:8*h; (0:h:8*h)*0];

[val,grad,hess] = helm1d.gshelm(rts,ejs,src,targ);

d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;
d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560	]/h^2;

err = d1*val - grad(5);
err2 = d2*val - hess(5);

%% part 3 - nonlocal helmholtz kernels (G_\phi)


beta = 1+0.5i;
gamma = 2;
[rts,ejs] = helm1d.find_roots(beta,gamma);

src = [];
src.r = [5;0];

h = 0.01;

targ = [];
targ.r = [0:h:8*h; (0:h:8*h)*0];

[val,grad,hess] = helm1d.gphihelm(rts,ejs,src,targ);

d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;
d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560	]/h^2;

err = d1*val - grad(5);
err2 = d2*val - hess(5);

