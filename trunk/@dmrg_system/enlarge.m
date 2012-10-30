function r = enlarge(varargin)
% ENLARGE - insert two sites in a dmrg system
%  ENLARGE(s,h,M,Ne) : enlarge system s using n.n. interaction h,
%    block dimensions M targeting Ne states
%  ENLARGE(s,h,M,Ne,opteigs) : enlarge system s using n.n. interaction h,
%    block dimensions M targeting Ne states using opteigs options for eigs.
%  ENLARGE(s,h,M,Ne,opteigs,B,t) : enlarge system s using
%    n.n. interaction h, block dimensions M, targeting Ne states using
%    opteigs options for eigs, and adding an external time dependent
%    magnetic field B at time t
%  INPUT:
%   s:  dmrg system
%   h:  n.n interaction between the two new sites (function handle)
%   M:  vector containing the dimensions of different blocks
%   B:  single-site time dep. magnetic field (function handle).
%  OUTPUT:
%   r:      dmrg system
%  NOTE: we do not test if the input arguments are of the correct type or
%   size!
  
  s = varargin{1};
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
  [Phi, E] = compute_phi(s,M,Ne,opteigs);
  [rhoL,rhoR] = PartialTrace2(Phi*Phi', M(s.l)*d, d*M(s.L-s.l-2));
  [VL,FL] = eig(rhoL);
  [VR,FR] = eig(rhoR);
  [FL,p] = sort(-diag(FL)); % order
  VL = VL(:,p); % reorder eigenvectors accordingly
  WL = VL(:,1:M(s.l+1)); % select first M(s.l+1)
  [FR,p] = sort(-diag(FR)); % order
  VR = VR(:,p); % reorder eigenvectors accordingly
  WR = VR(:,1:M(s.L-s.l-1)); % select first M(s.L-s.l-1)
  WLM = kron(WL,speye(d));
  WMR = kron(speye(d),WR);
  % insert 2 new sites in the middle and update hamiltonians
  % following definitions give correct hamiltonian for homogeneous h
  % and B only
  hLM = WLM'*kron(speye(M(s.l)),feval(h,s.l+1,s.l+2))*WLM;
  hMR = WMR'*kron(feval(h,s.l+3,s.l+4),speye(M(s.L-s.l-2)))*WMR;
  if nargin ~= 7
    HL = WL'*(kron(s.HL,speye(d))+s.hLM)*WL;
    hM = feval(h,s.l+2,s.l+3);
    HR = WR'*(s.hMR+kron(speye(d),s.HR))*WR;
  else
    HL = WL'*(kron(s.HL,speye(d))+s.hLM+kron(speye(M(s.l)),feval(B,t,s.l+1)))*WL;
    hM = feval(h,s.l+2,s.l+3)+kron(feval(B,t,s.l+2),speye(d))+...
         kron(speye(d),feval(B,t,s.l+3)); 
    HR = WR'*(kron(feval(B,t,s.l+4),speye(M(s.L-s.l-2)))+s.hMR+...
              kron(speye(d),s.HR))*WR; 
  end
  r = dmrg_system(s.L+2,s.l+1,HL,hLM,hM,hMR,HR,WL,WR);
  
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

