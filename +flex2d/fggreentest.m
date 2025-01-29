% testing greens function derivatives using finite differences 


src = [0; 0];
h = 0.05;
center = [1;1];
targs = [center(1)-4*h:h:center(1)+4*h; center(2) + zeros(1,9)];
rhoj = sqrt(2);

[rts,ejs] = flex2d.find_roots(3,-1);
[val,grad,hess,third,fourth] = flex2d.fggreen(src,targs,rts,ejs,true);

targs = [center(1) + zeros(1,9); center(2)-4*h:h:center(2)+4*h];
[valy,grady,hessy,thirdy,fourthy] = flex2d.fggreen(src,targs,rts,ejs,true);

dd = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280];

err1 = dd*val/h - grad(5,1,1);
err2 = dd*valy/h - grady(5,1,1);

err3 = dd*grad(:,:,1)/h - hess(5,1,1);
err4 = dd*grad(:,:,2)/h - hess(5,1,2);
err5 = dd*grady(:,:,2)/h - hess(5,1,3);

err6 = dd*hess(:,:,1)/h - third(5,1,1);
err7 = dd*hess(:,:,2)/h - third(5,1,2);
err8 = dd*hess(:,:,3)/h - third(5,1,3);
err9 = dd*hessy(:,:,3)/h - third(5,1,4);

err10 = dd*third(:,:,1)/h - fourth(5,1,1);
err11 = dd*third(:,:,2)/h - fourth(5,1,2);
err12 = dd*third(:,:,3)/h - fourth(5,1,3);
err13 = dd*third(:,:,4)/h - fourth(5,1,4);
err14 = dd*thirdy(:,:,4)/h - fourth(5,1,5);