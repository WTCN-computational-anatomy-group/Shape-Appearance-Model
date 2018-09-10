function varargout = Shoot(v0,settings,args)
% Geodesic shooting
% FORMAT psi = Shoot(v0,settings,args)
%
% v0       - Initial velocity field n1*n2*n3*3 (single prec. float)xi
% settings - 8 settings
%            - [1][2][3] Voxel sizes
%            - [4][5][6][7][8] Regularisation settings.
%            Regularisation uses the sum of
%            - [4] - absolute displacements
%            - [5] - laplacian
%            - [6] - bending energy
%            - [7] - linear elasticity mu
%            - [8] - linear elasticity lambda
% args     - Integration parameters
%            - [1] Num time steps
%
% psi      - Inverse deformation field n1*n2*n3*3 (single prec. float)
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
% initial momentum (u_0), computed (using differential operator L) by:
%     u_0 = L v_0
% Then (Ad_{\phi_t})^* m_0 is computed:
%     u_t = |d \phi_t| (d\phi_t)^T u_0(\phi_t)
% The velocity field at this time point is then obtained by using
% multigrid to solve:
%     v_t = L^{-1} u_t
%
% These equations can be found in:
% Younes (2007). "Jacobi fields in groups of diffeomorphisms and
% applications". Quarterly of Applied Mathematics, vol LXV,
% number 1, pages 113-134 (2007).
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
if nargin<3
    args = args0;
else
    if numel(args)<numel(args0)
        args = [args args0((numel(args)+1):end)];
    end
end

T     = args(1);   % # Time steps
d     = size(v0);
d     = d(1:3);
id    = Identity(d);

%if sum(v0(:).^2)==0
%    varargout{1} = id;
%    varargout{2} = v0;
%end

if ~isfinite(T)
    % Number of time steps from an educated guess about how far to move
    T = double(floor(sqrt(max(max(max(v0(:,:,:,1).^2+v0(:,:,:,2).^2+v0(:,:,:,3).^2)))))+1);
end

spm_diffeo('boundary',0);
F   = spm_shoot_greens('kernel',d,settings); % Could save time if this went outside
v   = v0;
u   = spm_diffeo('vel2mom',v,settings); % Initial momentum (u_0 = L v_0)
psi = id - v/T;

for t=2:abs(T)
    % The update of u_t is not exactly as described in the paper, but describing this might be a bit
    % tricky. The approach here was the most stable one I could find - although it does lose some
    % energy as < v_t, u_t> decreases over time steps.
    Jdp         = spm_diffeo('jacobian',id-v/T);
    u1          = zeros(size(u),'single');
    u1(:,:,:,1) = Jdp(:,:,:,1,1).*u(:,:,:,1) + Jdp(:,:,:,2,1).*u(:,:,:,2) + Jdp(:,:,:,3,1).*u(:,:,:,3);
    u1(:,:,:,2) = Jdp(:,:,:,1,2).*u(:,:,:,1) + Jdp(:,:,:,2,2).*u(:,:,:,2) + Jdp(:,:,:,3,2).*u(:,:,:,3);
    u1(:,:,:,3) = Jdp(:,:,:,1,3).*u(:,:,:,1) + Jdp(:,:,:,2,3).*u(:,:,:,2) + Jdp(:,:,:,3,3).*u(:,:,:,3);
    u           = spm_diffeo('pushc',u1,id+v/T);

    % v_t \gets L^g u_t
    v            = spm_shoot_greens(u,F,settings); % Convolve with Greens function of L

    % $\psi \gets \psi \circ (id - \tfrac{1}{T} v)$
    % Done in a slightly complicated way to make it easier to deal with wrapped boundary conditions
    % I found that simply using $\psi \gets \psi - \tfrac{1}{T} (D \psi) v$ was not so stable.
    dp           = id - v/T;
    psi(:,:,:,1) = spm_diffeo('bsplins',psi(:,:,:,1)-id(:,:,:,1),dp,[1 1 1  1 1 1]) + dp(:,:,:,1);
    psi(:,:,:,2) = spm_diffeo('bsplins',psi(:,:,:,2)-id(:,:,:,2),dp,[1 1 1  1 1 1]) + dp(:,:,:,2);
    psi(:,:,:,3) = spm_diffeo('bsplins',psi(:,:,:,3)-id(:,:,:,3),dp,[1 1 1  1 1 1]) + dp(:,:,:,3);
end
varargout{1} = psi;
varargout{2} = v;

