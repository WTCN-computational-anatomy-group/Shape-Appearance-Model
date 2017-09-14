function [T,iT,A] = OrthogonalisationMat(ZZ,S,WW,N,s)
Ka = s.Ka;
Kv = s.Kv;
K  = size(ZZ,1);
if Ka<=K || Kv<=K,
    ind   = {1:K};
else
    ind  = {1:Ka,(1:Kv)+Ka};
end

T   = zeros(K,K);
iT  = zeros(K,K);
for i=1:numel(ind)
    [T(ind{i},ind{i}),iT(ind{i},ind{i})] = OrthogonalisationMatrix(ZZ(ind{i},ind{i}),WW(ind{i},ind{i}),N);
end

ZZ1  = T*ZZ*T';
EZZ1 = T*(ZZ+S)*T';
WW1  = iT'*WW*iT;
q    = zeros(K,1)-0.5*log(N);
%q   = -log(N./diag(WW1))/2;
q    = min(max(q,-10),10);
Q    = diag(exp(q));
A    = SetReg(Q*EZZ1*Q,N,s);
E    = 0.5*(trace(Q*ZZ1*Q*A) + trace(WW1*inv(Q*Q)));
%fprintf('\n%d %g %g %g\n', 0, 0.5*trace(Q*ZZ1*Q*A), 0.5*trace(WW1*inv(Q*Q)), E)

for iter=1:100
    A   = SetReg(Q*EZZ1*Q,N,s);
    oE0 = E;

    for subit=1:10
        R  = A.*ZZ1'+A'.*ZZ1;
        g1 = Q*R*diag(Q);
        g2 =-2*(Q^2\diag(WW1));
        g  = g1+g2;

        H1 = Q*R*Q + diag(g1);
        H2 = 4*(Q^2\WW1);
        H  = H1+H2;

        H  = H + max(diag(H))*1e-8*eye(size(H));
        q  = q - H\g;
        q  = min(max(q,-10),10);
        Q  = diag(exp(q));

        oE = E;
        E  = 0.5*(trace(Q*ZZ1*Q*A) + trace(WW1*inv(Q*Q)));
        %fprintf('%d %g %g %g\n', iter, 0.5*trace(Q*ZZ1*Q*A), 0.5*trace(WW1*inv(Q*Q)), E)
        if (oE-E)/E < 1e-8, break; end
    end
    if abs(oE0-E)/E < 1e-7, break; end
end
T  = Q*T;
iT = iT/Q;



function [T,iT] = OrthogonalisationMatrix(ZZ,WW,N)
% Use SVD to orthognalise
% FORMAT [T,iT] = OrthogonalisationMatrix(ZZ,WW,N)
%
% ZZ - Z*Z' + S
% WW  - W'*L*W
% N   - Number of observations
%
% T   - Orthogonalisation matrix
% iT  - Inverse of T - roughly
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

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
iDz      = inv(Dz+max(abs(Dz(:)))*1e-9*eye(size(Dz)));
iDw      = inv(Dw+max(abs(Dw(:)))*1e-9*eye(size(Dw)));
T        = D*V'*iDz*Vz';
iT       = Vw*iDw*U;


% %Code for working out the gradients and Hessians
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



