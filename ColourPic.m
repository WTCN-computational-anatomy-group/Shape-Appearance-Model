function pic = ColourPic(a,likelihood)
% Construct a colour image for display
% FORMAT pic = ColourPic(a,likelihood)
%
% a          - Image data to work from
% likelihood - Either 'normal', 'binomial' or 'multinomial'.
% pic        - The resulting colour picture
%
% This function is used for debugging purposes
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$
d   = [size(a) 1 1];
s   = round((d(3)+1)/2);
if d(3)==121, s = 50; end
if d(3)==91,  s = 37; end
pic = a(:,:,s,:);
if nargin<2, likelihood = 'normal'; end
switch lower(likelihood)
case {'binomial','binary'}
    pic = 1./(1+exp(-pic));
case {'multinomial','categorical'}
    pic = SoftMax(pic);
end
pic = reshape(pic,[d(1:2) d(4)]);
if d(4)==1
    pic = repmat(pic,[1 1 3]);
elseif d(4)>3
    pic = pic(:,:,1:3);
else
    pic = cat(3,pic,zeros([d(1:2) (3-d(4))],'single'));
end
pic = uint8((round(min(max(rot90(pic),0),1)*255)));

