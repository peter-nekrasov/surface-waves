srcs = [0; 0];
targs = [1:0.05:1.45; 1+0*1:0.05:1.45];
rhoj = sqrt(2);

[val,grad,hess,third,fourth] = struveKdiffgreen(rhoj,targs,targs);