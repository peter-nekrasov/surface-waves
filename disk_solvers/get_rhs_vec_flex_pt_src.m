function v = get_rhs_vec_flex_pt_src(k,nu,n,R,src)

if nargin < 5
    src = [0; 0];
end

targ = [R; 0];

[u,up,upp,uppp] = hkdiffgreen(k,src,targ);
up = up(1,1,1);
upp = upp(1,1,1);
uppp = uppp(1,1,1);


v = [upp + nu/R*up - nu/R^2*n^2*u; 
    uppp + 1/R*upp - (1+n^2*(2-nu))/R^2*up + (3-nu)/R^3*n^2*u];

end