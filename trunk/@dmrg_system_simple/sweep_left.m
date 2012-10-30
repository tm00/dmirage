function r = sweep_left(varargin)
% SWEEP_LEFT - Perform one sweep step to the left
%  SWEEP_LEFT(s1,s2,Phi,M) : sweep system s1 to the left using
%    the left block from system s2, and the previously computed target
%    states Phi and block dimensions M
%  INPUT:
%   s1,s2: simple dmrg system
%   Phi:   state vectors
%   M:     vector containing the dimensions of different blocks
%  OUTPUT:
%   r:  simple dmrg system
%  NOTE: we do not test if the input arguments are of the correct type or
%   size!
  
  s1 = varargin{1};
  s2 = varargin{2};
  Phi = varargin{3};
  M = varargin{4};
  d = M(1);
  rhoR = PartialTrace(Phi*Phi',d*M(s1.L-2-s1.l));
  [VR,FR] = eig(rhoR);
  [FR,p] = sort(-diag(FR));
  VR = VR(:,p);
  WR = VR(:,1:M(s1.L-1-s1.l));
  r = dmrg_system_simple(s1.L,s1.l-1,s2.WL,WR);
  
function A1 = PartialTrace(A,d)
% PARTIALTRACE computes the partial trace of A on the second factor, only.
% A is a matrix on \C^d2 \otimes \C^d.
% The output A1 is the partial trace of A over \C^d2.

d2 = size(A)/d;
Id2 = speye(d2);
A1 = zeros(d); 
for k=1:d2
  A1 = A1 + kron(Id2(k,:),speye(d)) * A * kron( Id2(:,k),speye(d));
end