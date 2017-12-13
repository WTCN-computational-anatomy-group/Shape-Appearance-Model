function varargout = Shoot(v0,prm,args)
% Geodesic shooting
% FORMAT phi = spm_shoot3di(v0,prm,args)
% v0     - Initial velocity field n1*n2*n3*3 (single prec. float)xi
% prm  - 8 settings
%        - [1][2][3] Voxel sizes
%        - [4][5][6][7][8] Regularisation settings.
%          Regularisation uses the sum of
%          - [4] - absolute displacements
%          - [5] - laplacian
%          - [6] - bending energy
%          - [7] - linear elasticity mu
%          - [8] - linear elasticity lambda
% args   - Integration parameters
%          - [1] Num time steps
%
% theta  - Inverse deformation field n1*n2*n3*3 (single prec. float)
%
% This code generates inverse deformations from
% initial velocity fields by gedesic shooting.  See the work of Miller,
% Younes and others.
%
% LDDMM (Beg et al) uses the following evolution equation:
%     d\phi/dt = v_t(\phi_t)
% where a variational procedure is used to find the stationary solution
% for the time varying velocity field.
% In principle though, once the initial velocity is known, then the
% velocity at subsequent time points can be computed.  This requires
% initial momentum (m_0), computed (using differential operator L) by:
%     m_0 = L v_0
% Then (Ad_{\phi_t})^* m_0 is computed:
%     m_t = |d \phi_t| (d\phi_t)^T m_0(\phi_t)
% The velocity field at this time point is then obtained by using
% multigrid to solve:
%     v_t = L^{-1} m_t
%
% These equations can be found in:
% Younes (2007). "Jacobi fields in groups of diffeomorphisms and
% applications". Quarterly of Applied Mathematics, vol LXV,
% number 1, pages 113-134 (2007).
%
% Multigrid is currently used to obtain v_t = L^{-1} m_t, but
% this could also be done by convolution with the Greens function
% N = L^{-1} (see e.g. Bro-Nielson).
%
%________________________________________________________
% (c) Wellcome Trust Centre for NeuroImaging (2013-2017)

% John Ashburner
% $Id$

if isempty(v0)
    varargout{1} = [];
    varargout{2} = [];
    return;
end

args0 = 8;
if nargin<3,
    args = args0;
else
    if numel(args)<numel(args0),
        args = [args args0((numel(args)+1):end)];
    end
end

N     = args(1);   % # Time steps
d     = size(v0);
d     = d(1:3);
id    = Identity(d);

%if sum(v0(:).^2)==0
%    varargout{1} = id;
%    varargout{2} = v0;
%end

if ~isfinite(N),
    % Number of time steps from an educated guess about how far to move
    N = double(floor(sqrt(max(max(max(v0(:,:,:,1).^2+v0(:,:,:,2).^2+v0(:,:,:,3).^2)))))+1);
end

spm_field('boundary',0);
F     = spm_shoot_greens('kernel',d,prm); % Could save time if this went outside
vt    = v0;
mt    = spm_diffeo('vel2mom',vt,prm); % Initial momentum (m_0 = L v_0)
theta = id - vt/N;

for t=2:abs(N)
    Jdp          = spm_diffeo('jacobian',id-vt/N);
    mt1          = zeros(size(mt),'single');
    mt1(:,:,:,1) = Jdp(:,:,:,1,1).*mt(:,:,:,1) + Jdp(:,:,:,2,1).*mt(:,:,:,2) + Jdp(:,:,:,3,1).*mt(:,:,:,3);
    mt1(:,:,:,2) = Jdp(:,:,:,1,2).*mt(:,:,:,1) + Jdp(:,:,:,2,2).*mt(:,:,:,2) + Jdp(:,:,:,3,2).*mt(:,:,:,3);
    mt1(:,:,:,3) = Jdp(:,:,:,1,3).*mt(:,:,:,1) + Jdp(:,:,:,2,3).*mt(:,:,:,2) + Jdp(:,:,:,3,3).*mt(:,:,:,3);
    mt           = spm_diffeo('pushc',mt1,id+vt/N);

    vt           = spm_shoot_greens(mt,F,prm);
    dp           = id - vt/N;

    theta(:,:,:,1) = spm_diffeo('bsplins',theta(:,:,:,1)-id(:,:,:,1),dp,[1 1 1  1 1 1]) + dp(:,:,:,1);
    theta(:,:,:,2) = spm_diffeo('bsplins',theta(:,:,:,2)-id(:,:,:,2),dp,[1 1 1  1 1 1]) + dp(:,:,:,2);
    theta(:,:,:,3) = spm_diffeo('bsplins',theta(:,:,:,3)-id(:,:,:,3),dp,[1 1 1  1 1 1]) + dp(:,:,:,3);
end
varargout{1} = theta;
varargout{2} = vt;

