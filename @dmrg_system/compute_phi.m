function [Phi, E] = compute_phi(varargin)
% COMPUTE_PHI - Compute target states of dmrg system
%   compute_phi(s,M,Ne) : compute Ne target states and energies of system
%     s
%   compute_phi(s,M,Ne,opteigs) : compute Ne target states and energies
%     of system s with opteigs options for eigs

  s = varargin{1};
  M = varargin{2};
  Ne = varargin{3};
  switch nargin
   case 3
    opteigs = struct([]);
   case 4
    opteigs = varargin{4};
   otherwise
    error('Wrong number of input arguments');
  end
  d = M(1);

  H = @(v) Hv(v,s);
  [Phi,E] = eigs(H,M(s.l)*d^2*M(s.L-s.l-2),Ne,opteigs.sigma,opteigs);
  E = diag(E)';
  
function w =  Hv (v, s)
% HV matrix-vector multiplication with Hamiltonian
%  w = Hv(v, s) computes the matrix-vector multiplication of 
%    the Hamiltonian determined by the DMRG system s on the vector v
%    returning the vector w; 
%    v and w must be of dimension M(s.l)*N*N*M(s.L-s.l-2)
  
  
  % get the dimensions
  dL = size(s.HL,1);
  dLM = size(s.hLM,1);
  dM = size(s.hM,1);
  dMR = size(s.hMR,1);
  dR = size(s.HR,1);
  dTot = dL*dM*dR;
  % you'll need to play with tensor product indices to decipher this ...
  u1 = reshape(v, [dM*dR dL])*s.HL';
  w1 = reshape(u1, [dL*dM*dR 1]);
  
  u2 = reshape(v, [dMR dLM])*s.hLM';
  w2 = reshape(u2, [dLM*dMR 1]);
  
  u31 = reshape(v, [dR dL*dM])';
  u32 = s.hM*reshape(u31, [dM dL*dR]);
  u33 = reshape(u32, [dL*dM dR])';
  w3 = reshape(u33, [dL*dM*dR 1]);

  u4 = s.hMR*reshape(v, [dMR dLM]);
  w4 = reshape(u4, [dLM*dMR 1]);
  
  u5 = s.HR*reshape(v, [dR dL*dM]);
  w5 = reshape(u5, [dL*dM*dR 1]);
  
  w = w1 + w2 + w3 + w4 + w5;
