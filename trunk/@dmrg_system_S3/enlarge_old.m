function r = enlarge(varargin)
% ENLARGE - insert two sites in a dmrg system with conserved S3
%  ENLARGE(s,h,M,Ne) : enlarge system s using n.n. interaction h,
%    block dimensions M targeting Ne states
%  ENLARGE(s,h,M,Ne,opteigs) : enlarge system s using n.n. interaction h,
%    block dimensions M targeting Ne states using opteigs options for
%    eigs.
%  ENLARGE(s,h,M,Ne,opteigs,B,t) : enlarge system s using
%    n.n. interaction h, block dimensions M, targeting Ne states using
%    opteigs options for eigs, and adding an external time dependent
%    magnetic field B at time t
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
   case 7 
    opteigs = varargin{5};
    B = varargin{6};
    t = varargin{7};
   otherwise
    error('Wrong number of input arguments');
  end    
  d = M(1); % single-site dimension
  S3 = sparse(1:d,1:d,0.5*(d-1)-(0:d-1),d,d);
  [Phi, E] = compute_phi(s,M,Ne,opteigs);
  [rhoL,rhoR] = PartialTrace2(Phi*Phi', M(u.l)*d, d*M(u.L-u.l-2));
  S3L = kron(s.S3L,eye(d)) + kron(eye(M(u.l)),S3);
  S3R = kron(S3,eye(M(u.L-u.l-2))) + kron(eye(d),s.S3R);
  [aL,sL] = eig(full(S3L));
  [aR,sR] = eig(full(S3R));
  % choose eigenvectors of S3L and S3R with largest weight in reduced
  % density matrices
  [w,p] = sort(-diag(aL'*rhoL*aL)); % order
  aL = aL(:,p); % reorder eigenvectors accordingly
  WL = aL(:,1:M(u.l+1)); % select first M(u.l+1)
  [w,p] = sort(-diag(aR'*rhoR*aR)); % order
  aR = aR(:,p); % reorder eigenvectors accordingly
  WR = aR(:,1:M(u.L-u.l-1)); % select first M(u.L-u.l-2)
  WLM = kron(WL,eye(d));
  WMR = kron(eye(d),WR);
  % insert 2 new sites in the middle and update hamiltonians
  % following definitions give correct hamiltonian for homogeneous h
  % and B only
  if nargin ~= 7
    HL = WL'*(kron(u.HL,eye(d))+u.hLM)*WL;
    hLM = WLM'*kron(eye(M(u.l)),feval(h,u.l+1,u.l+2))*WLM;
    hM = feval(h,u.l+2,u.l+3);
    hMR = WMR'*kron(feval(h,u.l+3,u.l+4),eye(M(u.L-u.l-2)))*WMR;
    HR = WR'*(u.hMR+kron(eye(d),u.HR))*WR;
  else
    HL = WL'*(kron(u.HL,eye(d))+u.hLM+kron(eye(M(u.l)),feval(B,t,u.l+1)))*WL;
    hLM = WLM'*kron(eye(M(u.l)),feval(h,u.l+1,u.l+2))*WLM;
    hM = feval(h,u.l+2,u.l+3)+kron(feval(B,t,u.l+2),eye(d))+...
         kron(eye(d),feval(B,t,u.l+3));
    hMR = WMR'*kron(feval(h,u.l+3,u.l+4),eye(M(u.L-u.l-2)))*WMR;
    HR = WR'*(kron(feval(B,t,u.l+4),eye(M(u.L-u.l-2)))+u.hMR+...
              kron(eye(d),u.HR))*WR;
  end
  S3L = WL'*S3L*WL;
  S3R = WR'*S3R*WR;
  r = dmrg_system_S3(u.L+2,u.l+1,HL,hLM,hM,hMR,HR,WL,WR,S3L,S3R);
  
function [A1,A2] = PartialTrace2(A, d1, d2)
% PARTIALTRACE2 computes the partial trace of A
% A is a matrix on \C^d1 \otimes \C^d2 
% The output A1 is the partial trace of A over the second space, A2
% is the partial trace of A over the first space.

Id1 = eye(d1);
Id2 = eye(d2);
A1 = zeros(d1);
A2 = zeros(d2);
for k=1:d2
  A1 = A1 + kron(Id1, Id2(k,:)) * A * kron(Id1, Id2(:,k));
end
for k=1:d1
  A2 = A2 + kron(Id1(k,:), Id2) * A * kron(Id1(:,k), Id2);
end