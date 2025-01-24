function v = get_rhs_vec_gaussian(s,nu,n,R,src)

if nargin < 5
    src = [0; 0];
end

rho = sqrt(src(1)^2 + src(2)^2);
theta = atan2(src(2),src(1));

gn = exp(-R^2/(2*s^2));
gnp = -R*gn/s^2;
gnpp = gn*R^2/s^4 - gn/s^2;
gnppp = - gn*R^3/s^6 + 3*gn*R/s^4;


v = -[gnpp + nu/R*gnp - nu/R^2*n^2*gn; 
    gnppp + 1/R*gnpp - (1+n^2*(2-nu))/R^2*gnp + (3-nu)/R^3*n^2*gn];

%v(1) = v(1)*exp(-1i*n*theta);
%v(2) = v(2)*exp(-1i*n*theta);

end