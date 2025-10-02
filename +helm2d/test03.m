%% test for derivatives of greens functions
% part 1 - check helmholtz kernels
% 

k = 3;

src = [];
src.r = [5;0];

h = 0.01;

[X,Y] = meshgrid(0:h:8*h); 
targ = [];
targ.r = [X(:) Y(:)].';

d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;
d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560]/h^2;

[val,grad,hess] = chnk.helm2d.green(k,src.r,targ.r);
val = reshape(val,size(X));
gradx = reshape(grad(:,:,1),size(X));
grady = reshape(grad(:,:,2),size(X));
hessxx = reshape(hess(:,:,1),size(X));
hessxy = reshape(hess(:,:,2),size(X));
hessyy = reshape(hess(:,:,3),size(X));

err1 = d1*val(5,:).' - gradx(5,5)
err2 = d2*val(5,:).' - hessxx(5,5)


%% part 2 - nonlocal helmholtz kernels (G_S)

alpha = 2;
beta = 2;
gamma = 3;
[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);

src = [];
src.r = [1;0];

h = 0.01;

[X,Y] = meshgrid(0:h:8*h); 
targ = [];
targ.r = [X(:) Y(:)].';

[val,grad,hess] = helm2d.gshelm(rts,ejs,src,targ);
val = reshape(val,size(X));
gradx = reshape(grad(:,:,1),size(X));
grady = reshape(grad(:,:,2),size(X));
hessxx = reshape(hess(:,:,1),size(X));
hessxy = reshape(hess(:,:,2),size(X));
hessyy = reshape(hess(:,:,3),size(X));

d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;
d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560]/h^2;

err1 = d1*val(5,:).' - gradx(5,5)
err2 = d2*val(5,:).' - hessxx(5,5)

err3 = d1*val(:,5) - grady(5,5)
err4 = d2*val(:,5) - hessyy(5,5)

err5 = sum((d1.'*d1).*val,'all') - hessxy(5,5)

%% part 3 - nonlocal helmholtz kernels (G_\phi) % ~~ in progress ~~
% implement phi derivatives and come back to this

alpha = 0.5;
beta = -2;
gamma = 3;
[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);

src = [];
src.r = [3;0];

h = 0.01;

[X,Y] = meshgrid(0:h:8*h); 
targ = [];
targ.r = [X(:) Y(:)].';

[val,grad,hess] = helm2d.gphihelm(rts,ejs,src,targ);
val = reshape(val,size(X));
gradx = reshape(grad(:,:,1),size(X));
grady = reshape(grad(:,:,2),size(X));
hessxx = reshape(hess(:,:,1),size(X));
hessxy = reshape(hess(:,:,2),size(X));
hessyy = reshape(hess(:,:,3),size(X));

d1 = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280]/h;
d2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560]/h^2;

err1 = d1*val(5,:).' - gradx(5,5)
err2 = d2*val(5,:).' - hessxx(5,5)

err3 = d1*val(:,5) - grady(5,5)
err4 = d2*val(:,5) - hessyy(5,5)

err5 = sum((d1.'*d1).*val,'all') - hessxy(5,5) 
