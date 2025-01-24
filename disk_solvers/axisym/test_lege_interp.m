[x,w,u] = lege.exps(16);

f = -(x-1/2).^2;

cs = u*f;

[pols,~] = lege.pols(-1,15);

v0 = pols.'*u*f;