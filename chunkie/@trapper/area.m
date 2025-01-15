function a = area(trap)
%AREA compute area of a closed, 2D trapper curve

assert(trap.dim == 2,'area only well-defined for 2d curves');

wts = trap.wts;
rnorm = trap.n;
a = sum(sum(wts.*sum(rnorm.*(trap.r),1)))/trap.dim;

end
