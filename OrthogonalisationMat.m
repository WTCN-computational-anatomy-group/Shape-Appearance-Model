function [T,iT,A] = OrthogonalisationMat(ZZ,S,WW,N,s)
% Orthogonalisation matrix
% FORMAT [T,iT,A] = OrthogonalisationMat(ZZ,S,WW,N,s)
%
% ZZ    - Z*Z', where Z are the latent variables
% S     - E[Z*Z'] = Z*Z' + S
% WW    - Wa'*La*Wa + Wv'*Lv*Wv
% N     - size(Z,2)
% s     - Settings.
%
% T     - Transform.
% iT    - Inverse transform.
% A     - Precision matrix for distribution of latent variables.
%
% Compute a suitable orthogonalisation matrix (T) and its inverse (iT).
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

Ka  = s.Ka;
Kv  = s.Kv;
K   = size(ZZ,1);
T   = zeros(K,K);
iT  = zeros(K,K);

% Figure out which blocks of ZZ should be orthogonalised.
if Ka==K && Kv==K
    ind   = {1:K};
else
    ind  = {1:Ka,(1:Kv)+Ka};
end

% Find matrix T (and iT - its inverse) that orthogonalises the blocks
for i=1:numel(ind)
    [T(ind{i},ind{i}),iT(ind{i},ind{i})] = OrthogonalisationMatrix(ZZ(ind{i},ind{i}),WW(ind{i},ind{i}));
end

% The remainder of the code figures out how to rescale the transform T, such that
% the ``energy'' is minimised.
ZZ1  = T*ZZ*T';      % Transformed Z*Z'
EZZ1 = T*(ZZ+S)*T';  % Transformed E[Z*Z']
WW1  = N*iT'*WW*iT;  % Transformed W'*L*W

% Model parameters (q) are logs of the rescaling factors, which ensures
% they remain positive.
q    = zeros(K,1);   % It would be -0.5*log(N) if B = eye(K);
Q    = diag(exp(q)); % Diagonal rescaling matrix
A    = SetReg(Q*EZZ1*Q,N,s); % E-step
E    = 0.5*(trace(Q*ZZ1*Q*A) + trace(WW1/(Q*Q))); % Initial objective function

for iter=1:100
    % E-step: Compute precision matrix (A) of latent variables 
    %         using current estimate of the transform (T).
    A   = SetReg(Q*EZZ1*Q,N,s);

    % Current objective function for outer loop 
    oE0 = E;

    % Find the rescaling of T that minimises the "energy" using
    % the current estimate of the A matrix. 
    for subit=1:10
        R  = A.*ZZ1'+A'.*ZZ1;

        % Gradient
        g1 = Q*R*diag(Q);
        g2 =-2*(Q^2\diag(WW1));
        g  = g1+g2;

        % Hessian
        H1 = Q*R*Q + diag(g1);
        H2 = 4*(Q^2\WW1);
        H  = H1+H2;

        % Gauss-Newton update (with regularisation)
        H  = H + max(diag(H))*1e-6*eye(size(H));
        q  = q - H\g;
        q  = min(max(q,-10),10);
        Q  = diag(exp(q));

        % Check objective function for convergence
        oE = E;
        E  = 0.5*(trace(Q*ZZ1*Q*A) + trace(WW1/(Q*Q)));
        if (oE-E)/E < 1e-8, break; end
    end
    % If outer-loop objective function is unchannged, then done.
    if abs(oE0-E)/E < 1e-7, break; end
end
T  = Q*T;
iT = iT/Q;

%==========================================================================
%
%==========================================================================
function [T,iT] = OrthogonalisationMatrix(ZZ,WW)
% Use SVD to orthognalise
% FORMAT [T,iT] = OrthogonalisationMatrix(ZZ,WW)
%
% ZZ - Z*Z' + S
% WW  - W'*L*W
%
% T   - Orthogonalisation matrix
% iT  - Inverse of T - roughly
%
% Finds a transform that make ZZ and WW orthogonal, but the resulting
% transform needs to be rescaled such that an "energy" term is minimised.
%__________________________________________________________________________

if isempty(ZZ) && isempty(WW)
    T  = [];
    iT = [];
    return;
end

[Vz,Dz2] = svd(double(ZZ));
[Vw,Dw2] = svd(double(WW));
Dz       = diag(sqrt(diag(Dz2)+eps));
Dw       = diag(sqrt(diag(Dw2)+eps));
[U,D,V]  = svd(Dw*Vw'*Vz*Dz');
Dz       = Dz+max(abs(Dz(:)))*1e-8*eye(size(Dz));
Dw       = Dw+max(abs(Dw(:)))*1e-8*eye(size(Dw));
T        = D*V'*(Dz\Vz');
iT       = Vw*(Dw\U);


% % Code for working out the gradients and Hessians of the
% % for the energy minimisation
% q   = sym('q',[3,1],'real');
% Q   = diag(exp(q));
% A   = sym('a',[3,3],'real');
% ZZ1 = sym('x',[3,3],'real');
% y   = sym('y',[3,1],'real');
% WW1 = diag(y);
% E   = trace(Q*ZZ1*Q*A) + trace(WW1*inv(Q*Q));
% 
% pretty(simplify(diff(E,sym('q1')),1000))
% pretty(simplify(diff(diff(E,sym('q1')),sym('q2')),1000))
% pretty(simplify(diff(diff(E,sym('q1')),sym('q1')),1000))
% 
% g1 =  Q*(A.*ZZ1'+A'.*ZZ1)*diag(Q);
% g2 = -Q^2\diag(WW1)*2;
% g  =  g1+g2;
% H1 =  Q*(A'.*ZZ1 + A.*ZZ1')*Q +diag(g1);
% H2 =  4*WW1*Q^(-2);
% H  =  H1+H2;
% 
% d1  = simplify(g(1)  -diff(E,sym('q1')),1000)
% d11 = simplify(H(1,1)-diff(diff(E,sym('q1')),sym('q1')),1000)
% d12 = simplify(H(1,2)-diff(diff(E,sym('q1')),sym('q2')),1000)



