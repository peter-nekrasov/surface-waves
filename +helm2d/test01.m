%% Finite difference test for Green's functions
% Checking that the Green's functions satisfy the following relation:
%       1/2(\Delta + \beta)G_S + \gamma G_\phi = 0
% away from the diagonal (x \neq y)

h = 0.05;
xs = 1:h:(1+8*h);
[X,Y] = meshgrid(xs);
alpha = 0.5;
beta = -1.5;
gamma = 2.5;
R = sqrt(X.^2 + Y.^2);

src = [];
src.r = [0;0];
targ = [];
targ.r = [X(:) Y(:)].';

[rts,ejs] = helm2d.find_roots(alpha,beta,gamma);
gs = helm2d.gshelm(rts,ejs,src,targ);
gphi = helm2d.gphihelm(rts,ejs,src,targ);

gs = reshape(gs,size(X));
gphi = reshape(gphi,size(X));

% finite difference stencils
d2 = zeros(9, 1);
d2(1) = -1/560;
d2(2) = 8/315;
d2(3) = -1/5;
d2(4) = 8/5;
d2(5) = -205/72;
d2(6) = 8/5;
d2(7) = -1/5;
d2(8) = 8/315;
d2(9) = -1/560;
d2 = d2 / h^2;

lap = zeros(9,9);
lap(:,5) = d2;
lap(5,:) = lap(5,:) + d2.';

err = abs(0.5*alpha*sum(lap.*gs,'all') + 0.5*beta*gs(5,5) + gamma*gphi(5,5)) / max(abs(gs(:))) 