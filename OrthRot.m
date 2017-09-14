function [Rw,Rz] = OrthRot(EZZ,WW,N,s)
% Use SVD to orthognalise
% FORMAT [Rw,Rz] = OrthRot(EZZ,WW,N,s)
%
% EZZ - Z*Z'
% WW  - W'*L*W
% N   - Number of observations
% s   - Settings
%
% Rw  - Rotation to apply to W
% Rz  - Rotation to apply to Z
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if isempty(EZZ) && isempty(WW)
    Rw = [];
    Rz = [];
    return;
end

[Vz,Sz2]   = svd(EZZ);
[Vw,Sw2]   = svd(WW);
Sz         = diag(sqrt(diag(Sz2)+eps));
Sw         = diag(sqrt(diag(Sw2)+eps));
[Um,Sm,Vm] = svd(Sw*Vw'*Vz*Sz');
Tw         = (Vw/Sw)*Um;
Tz         = (Vz/Sz)*Vm;
sm         = diag(Sm);
sw         = ones(size(sm));

for it=1:1000
    Rz = Tz*diag(sm./sw);
    Rw = Tw*diag(sw);

    [EA,B] = SetReg(Rz'*EZZ*Rz,N,s);

    fz  = diag(EA).*sm.^2;
    fw  = diag(B);
    osw = sw;
    sw  = real((fz./fw).^(1/4));
    E   = 0.5*sum(fz./(sw.^2)) + 0.5*sum(fw.*(sw.^2));
    if sqrt(sum(log(sw./osw).^2)) < 1e-9, break; end
end

Rz = Tz*diag(sm./sw);
Rw = Tw*diag(sw);

