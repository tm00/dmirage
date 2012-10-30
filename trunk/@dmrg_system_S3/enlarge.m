function r = enlarge(varargin)
% ENLARGE - insert two sites in a dmrg system with conserved S3
%  ENLARGE(s,h,M,Ne) : enlarge system s using n.n. interaction h,
%    block dimensions M targeting Ne states
%  ENLARGE(s,h,M,Ne,opteigs) : enlarge system s using n.n. interaction h,
%    block dimensions M targeting Ne states using opteigs options for
%    eigs.
%  INPUT:
%   s:  dmrg system
%   h:  n.n interaction between the two new sites (function handle)
%   M:  vector containing the dimensions of different blocks
%  OUTPUT:
%   r:        dmrg system
%  NOTE: we do not test if the input arguments are of the correct type or
%   size!
%  NOTE: enlarging is really only meaningful for translation invariant
%   system, i.e., homogeneous h and B
  
  
  s = varargin{1};
  u = s.dmrg_system;
  h = varargin{2};
  M = varargin{3};
  Ne = varargin{4};
  switch nargin
   case 4
    opteigs = struct([]);
   case 5
    opteigs = varargin{5};
   otherwise
    error('Wrong number of input arguments');
  end    
  d = M(1); % single-site dimension
  S3 = sparse(1:d,1:d,0.5*(d-1)-(0:d-1),d,d);
  S3L = kron(s.S3L,speye(d)) + kron(speye(M(u.l)),S3);
  S3R = kron(S3,speye(M(u.L-u.l-2))) + kron(speye(d),s.S3R);
  
  [Phi, E] = compute_phi(s,M,Ne,s.mag,opteigs);
  [rhoL,rhoR] = PartialTrace2(Phi*Phi', M(u.l)*d, d*M(u.L-u.l-2));
  [VL,FL] = jointEig(rhoL,S3L);
  [FL,p] = sort(-FL);
  VL = VL(:,p); % reorder eigenvectors accordingly
  WL = VL(:,1:M(u.l+1)); % select first M(s.l+1)
   
  [VR,FR] = jointEig(rhoR,S3R);
  [FR,p] = sort(-FR);
  VR = VR(:,p); % reorder eigenvectors accordingly
  WR = VR(:,1:M(u.L-u.l-1)); % select first M(s.L-s.l-1)
  
  WLM = kron(WL,speye(d));
  WMR = kron(speye(d),WR);
  % insert 2 new sites in the middle and update hamiltonians
  HL = WL'*(kron(u.HL,speye(d))+u.hLM)*WL;
  hLM = WLM'*kron(speye(M(u.l)),feval(h,u.l+1,u.l+2))*WLM;
  hM = feval(h,u.l+2,u.l+3);
  hMR = WMR'*kron(feval(h,u.l+3,u.l+4),speye(M(u.L-u.l-2)))*WMR;
  HR = WR'*(u.hMR+kron(speye(d),u.HR))*WR;
  S3L = WL'*S3L*WL;
  S3R = WR'*S3R*WR;
  r = dmrg_system_S3(u.L+2,u.l+1,HL,hLM,hM,hMR,HR,WL,WR,S3L,S3R,s.mag);
  
function [A1,A2] = PartialTrace2(A, d1, d2)
% PARTIALTRACE2 computes the partial trace of A
% A is a matrix on \C^d1 \otimes \C^d2 
% The output A1 is the partial trace of A over the second space, A2
% is the partial trace of A over the first space.

Id1 = speye(d1);
Id2 = speye(d2);
A1 = zeros(d1);
A2 = zeros(d2);
for k=1:d2
  A1 = A1 + kron(Id1, Id2(k,:)) * A * kron(Id1, Id2(:,k));
end
for k=1:d1
  A2 = A2 + kron(Id1(k,:), Id2) * A * kron(Id1(:,k), Id2);
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

