function r = sweep_left(varargin)
% SWEEP_LEFT - Perform one sweep step to the left
%  SWEEP_LEFT(s1,s2,Phi,h,M,Ne) : sweep system s1 to the left using
%    the left block from system s2, and the
%    previously computed target states Phi, n.n. interaction h, block
%    dimensions M, targeting Ne states.
%  SWEEP_LEFT(s1,s2,Phi,h,M,Ne,opteigs) : sweep system s1 to the left
%    using the left block from system s2, and the
%    previously computed target states Phi, n.n. interaction h, block
%    dimensions M, new left block Hamiltonian HL and transformation matrix
%    WL, targeting Ne states using opteigs options for eigs.
%  SWEEP_LEFT(s1,s2,,Phi,h,M,Ne,opteigs,B,t) : sweep system s1 to the
%    left using the left block from system s2, and the
%    previously computed target states Phi, n.n. interaction
%    h, block dimensions M, new left block Hamiltonian HL and
%    transformation matrix WL, targeting Ne states using opteigs options
%    for eigs, and adding an external time dependent magnetic field B at
%    time t 
%  INPUT:
%   s1,s2: dmrg system
%   Phi:   state vectors
%   h:     n.n interaction between (function handle)
%   M:     vector containing the dimensions of different blocks
%   B:  single-site time dep. magnetic field (function handle).
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
   case 9
    opteigs = varargin{7};
    B = varargin{8};
    t = varargin{9};
   otherwise
    error('Wrong number of input arguments');
  end    
  d = M(1);
  rhoR = PartialTrace(Phi*Phi',d*M(s1.L-2-s1.l));
%  [rhoL,rhoR] = PartialTrace2(Phi*Phi', M(s1.l)*d, d*M(s1.L-2-s1.l));
  [VR,FR] = eig(rhoR);
  [FR,p] = sort(-diag(FR));
  VR = VR(:,p);
  WR = VR(:,1:M(s1.L-1-s1.l));
  WMR = kron(speye(d),WR);
  hMR = WMR'*kron(feval(h,s1.l+1,s1.l+2),speye(M(s1.L-2-s1.l)))*WMR;
  if s1.l == 2
    hLM = feval(h,s1.l-1,s1.l);
  else
    WLM = kron(s2.WL,speye(d));
    hLM = WLM'*kron(speye(M(s1.l-2)),feval(h,s1.l-1,s1.l))*WLM;
  end
  if nargin ~= 9
    hM = feval(h,s1.l,s1.l+1);
    HR = WR'*(s1.hMR+kron(speye(d),s1.HR))*WR;
  else
    hM = feval(h,s1.l,s1.l+1)+kron(feval(B,t,s1.l),speye(d))+...
         kron(speye(d),feval(B,t,s1.l+1));
    HR = WR'*(kron(feval(B,t,s1.l+2),speye(M(s1.L-s1.l-2)))+s1.hMR+...
              kron(speye(d),s1.HR))*WR;
  end
  r = dmrg_system(s1.L,s1.l-1,s2.HL,hLM,hM,hMR,HR,s2.WL,WR);
  
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