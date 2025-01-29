function submat= kern(zk,srcinfo,targinfo,type,varargin)
%FLEX2D.KERN standard Modified biharmonic layer potential kernels in 2D
% 
% Syntax: submat = flex2d.kern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels based on G(x,y) = i/4 H_0^{(1)}(zk |x-y|)
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
% D'(x,y) = \nabla_{n_x} \nabla_{n_y} G(x,y)
%
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (d/d2) it doesn't need to
%                be provided. sprime requires tangent info in
%                targinfo.d
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'dprime', normal derivative of double layer D'
%                type == 'c', combined layer kernel coef(1) D + coef(2) S
%                type == 'stau', tangential derivative of single layer
%                type == 'all', returns all four layer potentials, 
%                       [coef(1,1)*D coef(1,2)*S; coef(2,1)*D' coef(2,2)*S']
%                type == 'c2trans' returns the combined field, and the 
%                          normal derivative of the combined field
%                        [coef(1)*D + coef(2)*S; coef(1)*D' + coef(2)*S']
%                type == 'trans_rep' returns the potential corresponding
%                           to the transmission representation
%                        [coef(1)*D coef(2)*S]
%                type == 'trans_rep_prime' returns the normal derivative
%                          corresponding to the transmission representation
%                        [coef(1)*D' coef(2)*S']
%                type == 'trans_rep_grad' returns the gradient corresponding
%                         to the transmission representation
%                        [coef(1)*d_x D coef(2)*d_x S;
%                         coef(1)*d_y D coef(2)*d_y S]
%
%   varargin{1} - coef: length 2 array in the combined layer 
%                 formula, 2x2 matrix for all kernels
%                 otherwise does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also CHNK.HELM2D.GREEN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);



% clamped plate kernel for modified biharmonic problem
if strcmpi(type, 'clamped-plate')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   
   [~, ~, hess, third, ~] = flex2d.hkdiffgreen(zk, src, targ); 

   [~, ~, ~, ~, forth] = flex2d.hkdiffgreen(zk, src, targ, true);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;
   
   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   normaldot = nx.*nxtarg + ny.*nytarg;

   rn = rx.*nx + ry.*ny;
   rtau = rx.*taux + ry.*tauy;
   ntargtau = nxtarg.*taux + nytarg.*tauy;


   rntarg = rx.*nxtarg + ry.*nytarg;


   Kxx = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny))) - ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

   Kxy = -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

   Kyx = -(1/(2*zk^2).*(forth(:, :, 1).*(nx.*nx.*nx.*nxtarg) + forth(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
          forth(:, :, 3).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + forth(:, :, 4).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
          forth(:, :, 5).*(ny.*ny.*ny.*nytarg)) ) - ...
          (3/(2*zk^2).*(forth(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ forth(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
          forth(:, :, 3).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
          forth(:, :, 4).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) +...
          forth(:, :, 5).*(ny.*tauy.*tauy.*nytarg))) + ...
          1/pi.*(-3*rn.*rntarg./(r2.^2) + 4.*(rn.^3).*rntarg./(r2.^3) + 3*(rn.*rtau.*ntargtau)./ (r2.^2));
 %         1/pi.*(-3.*(rn.^2).*normaldot./(r2.^2) + 4.*(rn.^3).*rntarg./(r2.^3));


   Kyy = -(1/(2*zk^2).*(third(:,:, 1).*(nx.*nx.*nxtarg) +third(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + third(:, :, 3).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         third(:, :,4).*(ny.*ny.*nytarg))) + ...
         (1/(2*zk^2).*(third(:,:, 1).*(taux.*taux.*nxtarg) +third(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + third(:, :, 3).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         third(:, :,4).*(tauy.*tauy.*nytarg)));

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = Kxx;
  submat(1:2:end,2:2:end) = Kxy;
    
  submat(2:2:end,1:2:end) = Kyx;
  submat(2:2:end,2:2:end) = Kyy;
end



% clamped plate kernel for modified biharmonic problem
if strcmpi(type, 'clamped-plate K21')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   [~, ~, ~, ~, forth] = flex2d.bhgreen(src, targ);  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;
   
   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   normaldot = nx.*nxtarg + ny.*nytarg;

   rn = rx.*nx + ry.*ny;
   rtau = rx.*taux + ry.*tauy;
   ntargtau = nxtarg.*taux + nytarg.*tauy;

   rntarg = rx.*nxtarg + ry.*nytarg;

   submat = -((forth(:, :, 1).*(nx.*nx.*nx.*nxtarg) + forth(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
          forth(:, :, 3).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + forth(:, :, 4).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
          forth(:, :, 5).*(ny.*ny.*ny.*nytarg)) ) - ...
          (3*(forth(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ forth(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
          forth(:, :, 3).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
          forth(:, :, 4).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) +...
          forth(:, :, 5).*(ny.*tauy.*tauy.*nytarg))); 
end



% free plate kernel for the modified biharmonic problem (K11 with no
% hilbert transform subtraction, K12 kernel, K22 with no curvature part) K21 is
% handled in a separate type.
if strcmpi(type, 'free plate first part')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;
   coefs = varargin{1};

   [~, ~, hess, third, ~] = flex2d.hkdiffgreen(zk, src, targ);     
   [~, ~, ~, ~, forth] = flex2d.hkdiffgreen(zk, src, targ, true);   
   [~, ~, ~, ~, forthbh] = flex2d.bhgreen(src, targ);    
   % forth = forth + 2*zk^2*forthbh;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);


   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy); 

   taux = dx ./ ds;
   tauy = dy ./ ds;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   
   K11 = -(1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny))) - ...
       coefs(1)./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ;  % first kernel with no hilbert transforms (G_{nx nx ny + nu G_{taux taux ny}).

   K12 =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)+...
           coefs(1)/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}
   
   K21 = kappax./(2*zk^2).*(1-coefs(1)).*((third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
            third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
           (third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ) ...
        - 1/(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          forth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) - ...
          ((2-coefs(1))/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) ) ; % - ...          
          % (1+coefs(1))/(4*pi).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2));

   K22 = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))  +...
        (2-coefs(1))/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) + ... % G_{nx nx nx} + (2-nu) G_{taux taux nx}
        + kappax.*(1-coefs(1)).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg));

          
    
  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

end

% free plate kernel for the modified biharmonic problem (K11 with no
% hilbert transform subtraction, K12 kernel, K22 with no curvature part) K21 is
% handled in a separate type.
if strcmpi(type, 'free plate first part bh')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;
   coefs = varargin{1};

   % [~, ~, hess, third, ~] = flex2d.hkdiffgreen(zk, src, targ);     
   % [~, ~, ~, ~, forth] = flex2d.hkdiffgreen(zk, src, targ, true);   
   [~, ~, ~, ~, forthbh] = flex2d.bhgreen(src, targ);    
   forth = 2*zk^2*forthbh;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);


   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy); 

   taux = dx ./ ds;
   tauy = dy ./ ds;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   
   K11 = 0;

   K12 =  0;

   K21 = - 1/(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          forth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) - ...
          ((2-coefs(1))/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) ) - ...          
          (1+coefs(1))/(4*pi).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2));

   K22 = 0;
    
  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

end

% free plate kernel K21 first part (biharmonic part only) 
if strcmpi(type, 'free plate K21 first part bh')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

   [~, ~,~, ~, forth] = flex2d.bhgreen(src, targ);  % biharmonic part

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);


   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   tauxtargnsrc = tauxtarg.*nx + tauytarg.*ny;
   tauxtargntarg = tauxtarg.*nxtarg + tauytarg.*nytarg;


   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   normaldot = nx.*nxtarg + ny.*nytarg;

   rn = rx.*nx + ry.*ny;

   rntarg = rx.*nxtarg + ry.*nytarg;

   
   rtautarg = rx.*tauxtarg + ry.*tauytarg;


   submat  = -(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          forth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny))  - ...
          (2-coefs(1)).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny))  - ...           
          (1+coefs(1))/(4*pi).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2));
   % third kernel with singularity subtraction and H[delta']
   % 
 
end

% kernels in K11 with hilbert transform subtractions. 
% (i.e. beta*(G_{nx nx tauy} + 1/4 H + nu*(G_{taux taux tauy} + 1/4 H))

if strcmpi(type, 'free plate hilbert subtract')                                 
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

   [~,~,~, third, ~] = flex2d.hkdiffgreen(zk, src, targ, true);            % Hankel part

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   [~,grad] = chnk.lap2d.green(src,targ,true); 

   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);
  
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    


   submat = -((1+ coefs(1))/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*tauy))) + ((1+ coefs(1))/2).*0.25*hilb -...
       ((1+ coefs(1))/2)*coefs(1).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*tauy)));   % hilbert subtraction 

end




% kernels in K21 that are coupled with Hilbert transforms. 
if strcmpi(type, 'free plate hilbert')                                  
   srctang = srcinfo.d;
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;
   coefs = varargin{1};
    
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
  
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   [~,grad] = chnk.lap2d.green(src,targ,true); 
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

   [~, ~,~, third, ~] = flex2d.hkdiffgreen(zk, src, targ, true);            % Hankel part

   % [~, ~,~, thirdbh, ~] = flex2d.bhgreen(src, targ);            % Hankel part
   % third = third + 2*zk^2*thirdbh;

   [~, ~,~, ~, forth] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part

   K11 =  -(1+ coefs(1))/2*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) )  - ...
       (1+ coefs(1))/2*coefs(1).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy))) ;

   K21 = kappax.*(1-coefs(1)).*(-((1+ coefs(1))/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ coefs(1))/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)))) ...
        -(1+ coefs(1))/2.*(1/(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*taux) + ...
          forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*tauy + 3*nxtarg.*nxtarg.*nytarg.*taux) + ...
          forth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*tauy + 3*nxtarg.*nytarg.*nytarg.*taux) +...
          forth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*tauy + nytarg.*nytarg.*nytarg.*taux) +...
          forth(:, :, 5).*(nytarg.*nytarg.*nytarg.*tauy)) ) - ...
          ((2-coefs(1))/2)*(1+coefs(1)).*(1/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*taux) + ...
          forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*tauy + tauxtarg.*tauxtarg.*nytarg.*taux + 2*tauxtarg.*tauytarg.*nxtarg.*taux) + ...
          forth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*tauy + 2*tauxtarg.*tauytarg.*nxtarg.*tauy + tauytarg.*tauytarg.*nxtarg.*taux + 2*tauxtarg.*tauytarg.*nytarg.*taux) +...
          forth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*tauy + 2*tauxtarg.*tauytarg.*nytarg.*tauy + tauytarg.*tauytarg.*nytarg.*taux) +...
         forth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*tauy)) ) ;

  K12 = 0;

  K22 = 0;

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;
    
end


% kernels in K21 that are coupled with Hilbert transforms. 
if strcmpi(type, 'free plate hilbert bh')                                  
   srctang = srcinfo.d;
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;
   coefs = varargin{1};
    
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
  
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappax = numer ./ denom; % target curvature

   [~,grad] = chnk.lap2d.green(src,targ,true); 
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

   [~, ~,~, thirdbh, ~] = flex2d.bhgreen(src, targ);            % Hankel part
   third = 2*zk^2*thirdbh;

   K11 =  -(1+ coefs(1))/2*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) )  - ...
       (1+ coefs(1))/2*coefs(1).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + (1+ coefs(1))/2*(1+coefs(1)).*0.25*hilb;

   K21 = kappax.*(1-coefs(1)).*(-((1+ coefs(1))/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ coefs(1))/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy))))  ;

  K12 = 0;

  K22 = 0;

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;
    
end



% Updated part in K21 that is not coupled with Hilbert transform. 
%(i.e. (1-nu)*kappa*(-G_{nx nx ny} + G_{taux taux ny}). 

if strcmpi(type, 'free plate K21 second part')                                                         
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

    
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);


   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


  
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 


   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
   
   zkimag = (1i)*zk;

   [~,~,~, third, ~] = flex2d.helmdiffgreen(zk, src, targ);            % Hankel part
   [~,~,~, thirdK, ~] = flex2d.helmdiffgreen(zkimag, src, targ);     % modified bessel K part

   submat = ((1-coefs(1))/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
        (1-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nx) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*ny))) - ...
       ((1-coefs(1))/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) - ...
        (1-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*ny)));                                 % (1-nu)*(- G_{nx nx ny} + G_{ny taux taux})  
end

% Updated part in K21 that is coupled with Hilbert transform. 
%(i.e. (1-nu)*(beta*(G_{taux taux tauy} + 1/4 H) - beta*(G_{nx nx tauy} + 1/4 H)). 

if strcmpi(type, 'free plate K21 hilbert part bh')                                                         
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

   [~, ~, ~,third, ~] = flex2d.bhgreen(zk, src, targ);            % Hankel part
   third = 2*zk^2*third; 

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   [~,grad] = chnk.lap2d.green(src,targ,true); 
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    


   submat =  (1-coefs(1)).*(-((1+ coefs(1))/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ coefs(1))/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*tauy))))  ;  % + ((1+ coefs(1))/2).*0.25*hilb - ((1+ coefs(1))/2).*0.25*hilb)



end

if strcmpi(type, 'free plate K21 hilbert part')                                                         
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

   [~, ~, ~,third, ~] = flex2d.hkdiffgreen(zk, src, targ, true);            % Hankel part

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    

   submat =  (1-coefs(1)).*(-((1+ coefs(1))/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ coefs(1))/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*tauy))))  ;  % + ((1+ coefs(1))/2).*0.25*hilb - ((1+ coefs(1))/2).*0.25*hilb)



end



% Updated part in K22. (i.e. (1-nu)*(G_{taux taux} - G_{nx nx})

if strcmpi(type, 'free plate K22 second part')
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};



   [~, ~, hess, ~, ~] = flex2d.helmdiffgreen(zk, src, targ);            % Hankel part
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, hessK, ~, ~] = flex2d.helmdiffgreen(zkimag, src, targ);     % modified bessel K part
   
   

  

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


  
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 


   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    
   
   submat =  (1-coefs(1)).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hessK(:, :, 1).*tauxtarg.*tauxtarg + hessK(:, :, 2).*(2*tauxtarg.*tauytarg) + ...
           hessK(:, :, 3).*tauytarg.*tauytarg)-...
           (1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)-...
           1/(2*zk^2).*(hessK(:, :, 1).*nxtarg.*nxtarg + hessK(:, :, 2).*(2*nxtarg.*nytarg) + hessK(:, :, 3).*nytarg.*nytarg)));


end

if strcmpi(type, 'free plate hilbert unsubtract')                                 
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

   [~,~,~, third, ~] = flex2d.helmdiffgreen(zk, src, targ);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~,~, ~, thirdK, ~] = flex2d.helmdiffgreen(zkimag, src, targ);     % modified bessel K part

  
  

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
    


   submat = -((1+ coefs(1))/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*taux) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*tauy))) -...
       ((1+ coefs(1))/2)*coefs(1).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*tauy))); % + ((1+ coefs(1))/2)*coefs(1).*0.25*hilb;   % hilbert unsubtraction 

end


% free plate kernel for the interior modified biharmonic problem (K11 with no
% hilbert transform subtraction, K12 kernel, K22 with no curvature part) K21 is
% handled in a separate type.
if strcmpi(type, 'free plate first part interior')
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   



   [~, ~, hess, third, ~] = flex2d.helmdiffgreen(zk, src, targ);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, hessK, thirdK, ~] = flex2d.helmdiffgreen(zkimag, src, targ);     % modified bessel K part
   


   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   
   
   
   K11 = (1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
        1/(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nx) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*ny))) + ...
       coefs(1).*(1/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) - ...
        1/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*ny)));  % first kernel with no hilbert transforms (G_{nx nx ny + nu G_{taux taux ny}).
       
    
       

   K12 =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)-...
           1/(2*zk^2).*(hessK(:, :, 1).*nxtarg.*nxtarg + hessK(:, :, 2).*(2*nxtarg.*nytarg) + hessK(:, :, 3).*nytarg.*nytarg)+...
           coefs(1)/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           coefs(1)/(2*zk^2).*(hessK(:, :, 1).*tauxtarg.*tauxtarg + hessK(:, :, 2).*(2*tauxtarg.*tauytarg) + ...
           hessK(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}
   
   K21 = 0;




   K22 = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + thirdK(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg)+...
        thirdK(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + thirdK(:, :, 4).*(nytarg.*nytarg.*nytarg)) +...
        (2-coefs(1))/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) - ...
        (2-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + thirdK(:, :, 4).*(tauytarg.*tauytarg.*nytarg)); % G_{nx nx nx} + (2-nu) G_{taux taux nx}
    
  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

end

% free plate kernel K21 for the interior modified biharmonic problem. This part
% handles the singularity subtraction and swap the evaluation to its
% asymptotic expansions if the targets and sources are close. 

if strcmpi(type, 'free plate K21 first part interior')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

   [~, ~,~, ~, forth] = flex2d.helmdiffgreen(zk, src, targ, true);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, ~,~, forthK] = flex2d.helmdiffgreen(zkimag, src, targ, true);     % modified bessel K part 

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   tauxtargnsrc = tauxtarg.*nx + tauytarg.*ny;
   tauxtargntarg = tauxtarg.*nxtarg + tauytarg.*nytarg;


   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   normaldot = nx.*nxtarg + ny.*nytarg;

   rn = rx.*nx + ry.*ny;

   rntarg = rx.*nxtarg + ry.*nytarg;

   
   rtautarg = rx.*tauxtarg + ry.*tauytarg;


   submat  = (1/(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          forth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) -...
          1/(2*zk^2).*(forthK(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forthK(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forthK(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          forthK(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forthK(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny))) + ...
          ((2-coefs(1))/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) -...
          (2-coefs(1))/(2*zk^2).*(forthK(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forthK(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forthK(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forthK(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forthK(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny))) - ...
          (3/(2*pi).*(rn.*rntarg./(r2.^2)) + 3/(2*pi).*(rntarg.^2).*(normaldot)./(r2.^2) -...
            (2/pi).*(rn.*(rntarg.^3))./(r2.^3)) - ...
            (2-coefs(1)).*(-1/(2*pi).*(tauxtargnsrc).*(tauxtargntarg)./(r2) +...
            1/(pi).*(rn.*rtautarg.*tauxtargntarg)./(r2.^2) + ...
            1/(pi).*(rtautarg.*tauxtargnsrc.*rntarg)./(r2.^2) -....
            2/(pi).*(rn.*rntarg.*(rtautarg.*rtautarg))./(r2.^3) + ...
            1/(2*pi).*(rn.*rntarg)./(r2.^2)) +...           
          ((1+coefs(1))/(4*pi)).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2))+...
          (3/(4*pi)).*(nxtarg.*nx + nytarg.*ny)./r2  + (2-coefs(1))/(4*pi).*(nxtarg.*nx + nytarg.*ny)./r2 - ...
          (2-coefs(1))/(2*pi).*(nxtarg.*nx + nytarg.*ny).*(rx.*tauxtarg + ry.*tauytarg).*(rx.*tauxtarg + ry.*tauytarg)./(r2.^2);
   % third kernel with singularity subtraction and H[delta']
   % 
 
end

% Updated part in K21 (interior) that is not coupled with Hilbert transform. 
%(i.e. (1-nu)*(G_{nx nx ny} - G_{taux taux ny}). 

if strcmpi(type, 'free plate K21 second part interior')                                                         
   srcnorm = srcinfo.n;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};

    
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);


   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


  
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 


   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;
   
   zkimag = (1i)*zk;

   [~,~,~, third, ~] = flex2d.helmdiffgreen(zk, src, targ);            % Hankel part
   [~,~,~, thirdK, ~] = flex2d.helmdiffgreen(zkimag, src, targ);     % modified bessel K part

   submat = -((1-coefs(1))/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
        (1-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(nxtarg.*nxtarg.*nx) + thirdK(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        thirdK(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        thirdK(:, :, 4).*(nytarg.*nytarg.*ny))) + ...
       ((1-coefs(1))/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) - ...
        (1-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + thirdK(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        thirdK(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        thirdK(:, :, 4).*(tauytarg.*tauytarg.*ny)));                                 % (1-nu)*(- G_{nx nx ny} + G_{ny taux taux})  
end









% K11 kernel of the supported plate with 'Hilbert transform subtractions'.
if strcmpi(type, 'supported plate K11')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   



   [~, third, ~] = flex2d.helmdiffgreen(zk, src, targ, true);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, thirdK, ~] = flex2d.helmdiffgreen(zkimag, src, targ, true);     % modified bessel K part
   [~,grad] = chnk.lap2d.green(src,targ,true);                                 % lap kernels for the adjoint of Hilbert transform

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                               % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   
   stau = 2*(grad(:,:,1).*nytarg - grad(:,:,2).*nxtarg);                        % kernel of adjoint of Hilbert transform 
   
   c = (3-coefs(1))/(1 + coefs(1));

   submat = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*tauxtarg) + third(:, :, 2).*(nx.*nx.*tauytarg + 2*nx.*ny.*tauxtarg) +...
        third(:, :, 4).*(2*nx.*ny.*tauytarg + ny.*ny.*tauxtarg) +...
        third(:, :, 6).*(ny.*ny.*tauytarg)) - ...
        1/(2*zk^2).*(thirdK(:, :, 1).*(nx.*nx.*tauxtarg) + thirdK(:, :, 2).*(nx.*nx.*tauytarg + 2*nx.*ny.*tauxtarg) +...
        thirdK(:, :, 4).*(2*nx.*ny.*tauytarg + ny.*ny.*tauxtarg) +...
        thirdK(:, :, 6).*(ny.*ny.*tauytarg))) - (0.25*stau) + ...
        c.*(-(1/(2*zk^2).*(third(:, :, 1).*(taux.*taux.*tauxtarg) +...
        third(:, :, 2).*(taux.*taux.*tauytarg + 2*taux.*tauy.*tauxtarg) +...
        third(:, :, 4).*(2*taux.*tauy.*tauytarg + tauy.*tauy.*tauxtarg) +...
        third(:, :, 6).*(tauy.*tauy.*tauytarg)) - ...
        1/(2*zk^2).*(thirdK(:, :, 1).*(taux.*taux.*tauxtarg) +...
        thirdK(:, :, 2).*(taux.*taux.*tauytarg + 2*taux.*tauy.*tauxtarg) +...
        thirdK(:, :, 4).*(2*taux.*tauy.*tauytarg + tauy.*tauy.*tauxtarg) +...
        thirdK(:, :, 6).*(tauy.*tauy.*tauytarg))) - (0.25)*stau);      
   
   % (G_{ny ny nx} - 0.25*Stau) + c*(G_{tauy tauy taux}-0.25 Stau)).
       


end


if strcmpi(type, 'test kernel')                                                         % test K21 kernel
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};



   [~, ~,~, ~, forth] = flex2d.helmdiffgreen(zk, src, targ, true);            % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, ~, ~, forthK] = flex2d.helmdiffgreen(zkimag, src, targ, true);     % modified bessel K part
    

  
   

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);


   ds = sqrt(dx.*dx+dy.*dy);
   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;
   tauxtargnsrc = tauxtarg.*nx + tauytarg.*ny;
   tauxtargntarg = tauxtarg.*nxtarg + tauytarg.*nytarg;

   normaldot = nx.*nxtarg + ny.*nytarg;
   tangentdot = taux.*tauxtarg + tauy.*tauytarg;

   rn = rx.*nx + ry.*ny;

   rntarg = rx.*nxtarg + ry.*nytarg;

   rtau = rx.*taux + ry.*tauy;

   rtautarg = rx.*tauxtarg + ry.*tauytarg;
   %inds = and((dists <= R),(dists > 0));

   submat  = -(1/(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          forth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) -...
          1/(2*zk^2).*(forthK(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + forthK(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          forthK(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          forthK(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          forthK(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny))) - ...
          ((2-coefs(1))/(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) -...
          (2-coefs(1))/(2*zk^2).*(forthK(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + forthK(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          forthK(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          forthK(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          forthK(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny))) + ...
          3/(2*pi).*(rn.*rntarg./(r2.^2)) + 3/(2*pi).*(rntarg.^2).*(normaldot)./(r2.^2) -...
            (2/pi).*(rn.*(rntarg.^3))./(r2.^3) + ...
            (2-coefs(1)).*(-1/(2*pi).*(tauxtargnsrc).*(tauxtargntarg)./(r2) +...
            1/(pi).*(rn.*rtautarg.*tauxtargntarg)./(r2.^2) + ...
            1/(pi).*(rtautarg.*tauxtargnsrc.*rntarg)./(r2.^2) -....
            2/(pi).*(rn.*rntarg.*(rtautarg.*rtautarg))./(r2.^3) + ...
            1/(2*pi).*(rn.*rntarg)./(r2.^2)) -...           
          ((1+coefs(1))/(4*pi)).*((taux.*tauxtarg + tauy.*tauytarg)./(r2) - 2*(rx.*tauxtarg + ry.*tauytarg).*(rx.*taux + ry.*tauy)./(r2.^2))-...
          (3/(4*pi)).*(nxtarg.*nx + nytarg.*ny)./r2  - (2-coefs(1))/(4*pi).*(nxtarg.*nx + nytarg.*ny)./r2 + ...
          (2-coefs(1))/(2*pi).*(nxtarg.*nx + nytarg.*ny).*(rx.*tauxtarg + ry.*tauytarg).*(rx.*tauxtarg + ry.*tauytarg)./(r2.^2);
          %+...
          % 3/(2*pi).*(rn.*rntarg./(r2.^2)) + 3/(2*pi).*(rntarg.^2).*(normaldot)./(r2.^2) -...
          % (2/pi).*(rn.*(rntarg.^3))./(r2.^3) + ...
          % (2-coefs(1)).*(-1/(2*pi).*(tauxtargnsrc).*(tauxtargntarg)./(r2) +...
          %  1/(pi).*(rn.*rtautarg.*tauxtargntarg)./(r2.^2) + ...
          % 1/(pi).*(rtautarg.*tauxtargnsrc.*rntarg)./(r2.^2) -....
          % 2/(pi).*(rn.*rntarg.*(rtautarg.*rtautarg))./(r2.^3) + ...
          % 1/(2*pi).*(rn.*rntarg)./(r2.^2));
   % third kernel with singularity subtraction and H[delta']
    % +...
    %(2-coefs(1))/(4*pi).*(nxtarg.*nx + nytarg.*ny)./r2 - ...
   % ((2-coefs(1))/(2*pi)).*(nxtarg.*nx + nytarg.*ny).*(rx.*tauxtarg + ry.*tauytarg).*(rx.*tauxtarg + ry.*tauytarg)./(r2.^2)
     % if any(inds(:))                                                      % if statement for evalauation points that are closed to the boundary
   %     temp = 1i*imag(submat(inds));
   % 
   %     nx = nx(inds); ny = ny(inds); 
   % 
   %     nxtarg = nxtarg(inds);
   %     nytarg = nytarg(inds);
   % 
   %     dx = dx(inds);
   %     dy = dy(inds);
   % 
   %     dx1 = dx1(inds);
   %     dy1 = dy1(inds);
   % 
   %     ds = sqrt(dx.*dx+dy.*dy);
   %     ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
   % 
   %     taux = dx./ds;                                                                       % normalization
   %     tauy = dy./ds;
   % 
   %     tauxtarg = dx1./ds1;
   %     tauytarg = dy1./ds1;
   % 
   %     rx = rx(inds);
   %     ry = ry(inds);
   %     r2 = r2(inds);
   % 
   %     tauxtargnsrc = tauxtarg.*nx + tauytarg.*ny;
   %     tauxtargntarg = tauxtarg.*nxtarg + tauytarg.*nytarg;
   % 
   %     normaldot = nx.*nxtarg + ny.*nytarg;
   %     tangentdot = taux.*tauxtarg + tauy.*tauytarg;
   % 
   %     rn = rx.*nx + ry.*ny;
   % 
   %     rntarg = rx.*nxtarg + ry.*nytarg;
   % 
   %     rtau = rx.*taux + ry.*tauy;
   % 
   %     rtautarg = rx.*tauxtarg + ry.*tauytarg;
   % 
   % 
   % 
   %     [smoothfirst, smoothsecond] = chnk.helm2d.smooth_part_derivatives(zk, rx, ry);
   % 
   %     eulerconstantpart = (smoothfirst(:, :, 1)).*(nxtarg.*nxtarg.*nxtarg.*nx) +...
   %         (smoothfirst(:, :, 2)).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
   %        (smoothfirst(:, :, 3)).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
   %        (smoothfirst(:, :, 4)).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
   %        (smoothfirst(:, :, 5)).*(nytarg.*nytarg.*nytarg.*ny)+ ...
   %        (2-coefs(1)).*((smoothfirst(:, :, 1)).*(tauxtarg.*tauxtarg.*nxtarg.*nx) +...
   %        (smoothfirst(:, :, 2)).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
   %        (smoothfirst(:, :, 3)).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
   %        (smoothfirst(:, :, 4)).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
   %        (smoothfirst(:, :, 5)).*(tauytarg.*tauytarg.*nytarg.*ny));
   % 
   %     puresmoothpart = (smoothsecond(:, :, 1)).*(nxtarg.*nxtarg.*nxtarg.*nx) +...
   %         (smoothsecond(:, :, 2)).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
   %        (smoothsecond(:, :, 3)).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
   %        (smoothsecond(:, :, 4)).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
   %        (smoothsecond(:, :, 5)).*(nytarg.*nytarg.*nytarg.*ny)+ ...
   %        (2-coefs(1)).*((smoothsecond(:, :, 1)).*(tauxtarg.*tauxtarg.*nxtarg.*nx) +...
   %        (smoothsecond(:, :, 2)).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
   %        (smoothsecond(:, :, 3)).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
   %        (smoothsecond(:, :, 4)).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
   %        (smoothsecond(:, :, 5)).*(tauytarg.*tauytarg.*nytarg.*ny));
   % 
   %     submat(inds) = 3/(2*pi).*(rn.*rntarg./(r2.^2)) + 3/(2*pi).*(rntarg.^2).*(normaldot)./(r2.^2) -...
   %          (2/pi).*(rn.*(rntarg.^3))./(r2.^3) + ...
   %          (2-coefs(1)).*(-1/(2*pi).*(tauxtargnsrc).*(tauxtargntarg)./(r2) +...
   %          1/(pi).*(rn.*rtautarg.*tauxtargntarg)./(r2.^2) + ...
   %          1/(pi).*(rtautarg.*tauxtargnsrc.*rntarg)./(r2.^2) -....
   %          2/(pi).*(rn.*rntarg.*(rtautarg.*rtautarg))./(r2.^3) + ...
   %          1/(2*pi).*(rn.*rntarg)./(r2.^2)) - 3/(4*pi).*normaldot./r2 +...
   %          ((2-coefs(1))/(2*pi)).*((rtautarg.*rtautarg).*(normaldot))./(r2.^2) -...
   %          ((2-coefs(1))/(4*pi)).*(normaldot)./r2 - ...
   %        ((1+coefs(1))/(4*pi)).*((tangentdot)./(r2) - 2*(rtautarg).*(rtau)./(r2.^2)) + temp + eulerconstantpart + puresmoothpart;
   % end


end                                                         

if strcmpi(type, "test kernel jump")
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   

   [~, ~, ~, third, ~] = flex2d.helmdiffgreen(zk, src, targ);           % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   % nxtarg = repmat((targnorm(1,:)).',1,ns);
   % nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;
   [~, ~, ~, thirdK, ~] = flex2d.helmdiffgreen(zkimag, src, targ);     % modified bessel K part

   % [~,grad] = chnk.lap2d.green(src,targ,true);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;


   % dx1 = repmat((targtang(1,:)).',1,ns);
   % dy1 = repmat((targtang(2,:)).',1,ns);
   % 
   % 
   % 
   % ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   

   % tauxtarg = dx1./ds1;
   % tauytarg = dy1./ds1;

   submat = -(1./(2*zk^2).*(third(:, :, 1).*(nx.*nx.*taux) + third(:, :, 2).*(nx.*nx.*tauy+ 2*nx.*ny.*taux) +...
        third(:, :, 3).*(2*nx.*ny.*tauy +ny.*ny.*taux) +...
        third(:, :, 4).*(ny.*ny.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nx.*nx.*taux) + thirdK(:, :, 2).*(nx.*nx.*tauy+ 2*nx.*ny.*taux) +...
        thirdK(:, :, 3).*(2*nx.*ny.*tauy +ny.*ny.*taux) +...
        thirdK(:, :, 4).*(ny.*ny.*tauy))) +...
       (1./(2*zk^2).*(third(:, :, 1).*(taux.*taux.*taux) + third(:, :, 2).*(taux.*taux.*tauy + 2*taux.*tauy.*taux) +...
       third(:, :, 3).*(2*taux.*tauy.*tauy + tauy.*tauy.*taux) +...
        third(:, :, 4).*(tauy.*tauy.*tauy)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(taux.*taux.*taux) + thirdK(:, :, 2).*(taux.*taux.*tauy + 2*taux.*tauy.*taux) +...
       thirdK(:, :, 3).*(2*taux.*tauy.*tauy + tauy.*tauy.*taux) +...
        thirdK(:, :, 4).*(tauy.*tauy.*tauy)));
end   










% first-kernel of evaluation for the clamped plate
if strcmpi(type, 'first kernel')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;



    [~, ~, ~, third, ~] = flex2d.hkdiffgreen(zk, src, targ);           % Hankel part
    
    zkimag = 1i*zk;
    %[~,~, ~, thirdK, ~] = flex2d.helmdiffgreen(zkimag, src, targ);     % modified bessel K part



    submat = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) ) - ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

end


% second-kernel of evaluation for the clamped plate
if strcmpi(type, 'second kernel')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);

    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    [~, ~, hess] = flex2d.hkdiffgreen(zk, src, targ);           % Hankel part
    
    submat =  -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}
end


if strcmpi(type, 'free plate eval first')                                               % G_{ny}
   srcnorm = srcinfo.n;
   [~,grad] = flex2d.helmdiffgreen(zk,src,targ);        % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   


   zkimag = (1i)*zk;
   [~,gradK] = flex2d.helmdiffgreen(zkimag,src,targ);    % modified bessel K part

 
  

   submat = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny) + ...
       1/(2*zk^2).*(gradK(:, :, 1).*(nx) + gradK(:, :, 2).*ny)); 


end

if strcmpi(type, 'free plate eval first hilbert')                                               % G_{tauy} (supposed to coupled with the hilbert transform)
   srctang = srcinfo.d;
   coefs = varargin{1};

   [~,grad] = flex2d.helmdiffgreen(zk,src,targ);        % Hankel part

  
   zkimag = (1i)*zk;
   [~,gradK] = flex2d.helmdiffgreen(zkimag,src,targ);    % modified bessel K part

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);
 

   taux = dx./ds;                                                                       % normalization
   tauy = dy./ds;



    submat = ((1 + coefs(1))/2).*(-1/(2*zk^2).*(grad(:, :, 1).*(taux) + grad(:, :, 2).*tauy) + ...
       1/(2*zk^2).*(gradK(:, :, 1).*(taux) + gradK(:, :, 2).*tauy));                    % G_{tauy}
end


if strcmpi(type, 'free plate eval second')                                          % G = 1/(2k^2) (i/4 H_0^{1} - 1/2pi K_0)

   [val,~] = flex2d.helmdiffgreen(zk,src,targ);        % Hankel part

   zkimag = (1i)*zk;
   [valK,~] = flex2d.helmdiffgreen(zkimag,src,targ);    % modified bessel K part

   submat = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;

end



if strcmpi(type, 'fluid problem first evaluation')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   
   [~, ~, hess, ~, ~] = flex2d.helmdiffgreen(zk, src, targ);           % Hankel part
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   zkimag = (1i)*zk;
   [~, ~, hessK, ~, ~] = flex2d.helmdiffgreen(zkimag, src, targ);     % modified bessel K part

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;
   
   submat = (1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))-...
           1/(2*zk^2).*(hessK(:, :, 1).*(nx.*nx) + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*(ny.*ny)))-...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))-...
           1/(2*zk^2).*(hessK(:, :, 1).*(taux.*taux) + hessK(:, :, 2).*(2*taux.*tauy) + hessK(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}
end

if strcmpi(type, 'fluid problem second evaluation hilbert')
   srctang = srcinfo.d;
   
   [~, ~, hess, ~, ~] = flex2d.helmdiffgreen(zk, src, targ);           % Hankel part

   zkimag = (1i)*zk;
   [~, ~, hessK, ~, ~] = flex2d.helmdiffgreen(zkimag, src, targ);     % modified bessel K part

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   submat = (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))-...
           1/(2*zk^2).*(hessK(:, :, 1).*(taux.*taux) + hessK(:, :, 2).*(2*taux.*tauy) + hessK(:, :, 3).*(tauy.*tauy))); %  G_{tauy tauy}
end




end


