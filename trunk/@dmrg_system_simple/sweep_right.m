function r = sweep_right(varargin)
% SWEEP_RIGHT - Perform one sweep step to the right
%  SWEEP_RIGHT(s1,s2,Phi,M) : sweep system s1 to the right using
%    the right block from system s2, and the previously computed target
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
  rhoL = PartialTrace(Phi*Phi',d*M(s1.l));
  [VL,FL] = eig(rhoL);
  [FL,p] = sort(-diag(FL));
  VL = VL(:,p);
  WL = VL(:,1:M(s1.l+1));
  WLM = kron(WL,speye(d));
  r = dmrg_system_simple(s1.L,s1.l+1,WL,s2.WR);


function A1 = PartialTrace(A,d)
% PARTIALTRACE computes the partial trace of A on the first factor, only.
% A is a matrix on \C^d \otimes \C^d2.
% The output A1 is the partial trace of A over \C^d2.

d2 = size(A)/d;
Id2 = speye(d2);
A1 = zeros(d); 
for k=1:d2
  A1 = A1 + kron(speye(d), Id2(k,:)) * A * kron(speye(d), Id2(:,k));
end

