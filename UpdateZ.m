function [z,S,L,omisc] = UpdateZ(z,f,mu,Wa,Wv,A,s,noise)
% Update latent variables for one observation
% FORMAT [z,S,L,omisc] = UpdateZ(z,f,mu,Wa,Wv,A,s,noise)
%
% z     - Vector of latent variables
% f     - This observation
% mu    - Mean
% Wa    - Appearance basis functions
% Wv    - Shape basis functions
% A     - Precision matrix of z (assumed zero mean)
% s     - Settings. Uses s.v_settings, s.nit, s.vx & s.int_args
% noise - Noise precision (Gaussian model only)
%
% S     - Covariance of uncertainty of z (Laplace approximation)
% L     - Log-likelihood
% omisc - An assortment of statistics
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

[z,S,L,omisc] = UpdateZpar({z(:)},{f},mu,Wa,Wv,A,s,noise);
z = z{1};
S = S{1};
