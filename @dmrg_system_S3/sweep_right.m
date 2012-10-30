function r = sweep_right(varargin)
% SWEEP_RIGHT - Perform one sweep step to the right
%  SWEEP_RIGHT(s1,s2,Phi,h,M,Ne) : sweep system s1 to the right using
%    the right block from system s2, and the
%    previously computed target states Phi, n.n. interaction h, block
%    dimensions M, targeting Ne states.
%  SWEEP_RIGHT(s1,s2,Phi,h,M,Ne,opteigs) : sweep system s to the right
%    using the right block from system s2, and the
%    previously computed target states Phi, n.n. interaction h, block
%    dimensions M, new right block Hamiltonian HR and transformation matrix
%    WR, targeting Ne states using opteigs options for eigs.
%  INPUT:
%   s1,s2: dmrg system
%   Phi:   state vectors
%   h:     n.n interaction between (function m-file)
%   M:     vector containing the dimensions of different blocks
%  OUTPUT:
%   r:  dmrg system
%  NOTE: we do not test if the input arguments are of the correct type or
%   size!
  
  s1 = varargin{1};
  s2 = varargin{2};
  Phi = varargin{3};
  h = varargin{4};
  M = varargin{5};
  Ne = varargin{6};
  switch nargin
   case 6
    opteigs = struct([]);
   case 7
    opteigs = varargin{7};
   otherwise
    error('Wrong number of input arguments');
  end    
  d = M(1);
  u1 = s1.dmrg_system;
  u2 = s2.dmrg_system;
  S3 = sparse(1:d,1:d,0.5*(d-1)-(0:d-1),d,d);
  rhoL = PartialTrace(Phi*Phi',d*M(u1.l));
  S3L = kron(s1.S3L,speye(d)) + kron(speye(M(u1.l)),S3);
  [VL,FL] = jointEig(rhoL,S3L);
  [FL,p] = sort(-FL);
  VL = VL(:,p); % reorder eigenvectors accordingly
  WL = VL(:,1:M(u1.l+1)); % select first M(s.l+1)
  WLM = kron(WL,speye(d));
  hLM = WLM'*kron(speye(M(u1.l)),feval(h,u1.l+1,u1.l+2))*WLM;
  if u1.l == u1.L-4
    hMR = feval(h,u1.l+3,u1.l+4);
  else
    WMR = kron(speye(d),u2.WR);
    hMR = WMR'*kron(feval(h,u1.l+3,u1.l+4),speye(M(u1.L-u1.l-4)))*WMR;
  end
  hM = feval(h,u1.l+2,u1.l+3);
  HL = WL'*(kron(u1.HL,speye(d))+ u1.hLM)*WL;
  S3L = WL'*S3L*WL;
  r = dmrg_system_S3(u1.L,u1.l+1,HL,hLM,hM,hMR,u2.HR,WL,u2.WR,S3L,s2.S3R,s1.mag);
  
 
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

function [V,F] = jointEig(rho, S3)
% JOINTEIG eigenvalue decomposition of rho such that the eigenvectors are
% also eigenvectors of S3
  
% S3 should be a diagonal matrix of half integers and rho should commute
% with S3 (both upto finite precision)
[I,J] = find(S3);
if I~=J
  error('S3 is not diagonal!')
end
% round values in s to the nearest half-integer, everything else is
% finite precision errors
s = full(diag(S3));
s = round(s./0.5).*0.5;
[I,J,K] = find(rho>eps);
% [rho,S3] = 0 iff (s_i - s_j)rho_(ij) = 0 for all i,j
if (~isempty(find((s(I)-s(J)).*K > eps)))
  error('commutator of rho and S3 is not zero!')
end
% we are in business now
% first produce a list x of S3 values without multiplicities
x = [];
y = s;
while (length(y)~=0)
  x(end+1)=max(y);
  y=y(find(y~=max(y)));
end
% eigenvectors of S3 are given by standard basis
Id = eye(size(S3));
% initialize eigenvectors and eigenvalues of rho
V = zeros(size(rho));
F = zeros(length(rho),1);
% diagonalize rho within each subspace separately, note this solves
% several small eigenvalue problems instead of one big one
for k=1:length(x)
  ind = find(s == x(k));
  % restrict rho to subspace and diagonalize
  rhos = Id(:,ind)'*rho*Id(:,ind);
  [Vs,Fs] = eig(rhos);
  % put the result back into the large space
  V(:,ind) = Id(:,ind)*Vs;
  F(ind) = diag(Fs);
end