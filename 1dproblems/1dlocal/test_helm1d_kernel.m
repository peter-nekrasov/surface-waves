clear 
close all
addpath(genpath('..'))

src = [];
src.r = [5; 0];

targs = [];
targs.r = [0:0.02:10; (0:0.02:10)*0];

[val,grad] = helm1d.green(8,src,targs);

figure(1);
plot(targs.r(1,:),real(val),targs.r(1,:),imag(val));

figure(2);
plot(targs.r(1,:),real(grad),targs.r(1,:),imag(grad));

rmpath(genpath('..'))
