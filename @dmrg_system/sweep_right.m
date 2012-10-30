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
%  SWEEP_RIGHT(s1,s2,,Phi,h,M,Ne,opteigs,B,t) : sweep system s1 to the
%    right using the right block from system s2, and the
%    previously computed target states Phi, n.n. interaction h,
%    block dimensions M, new right block Hamiltonian HR and transformation
%    matrix WR, targeting Ne states using opteigs options for eigs,  and
%    adding an external time dependent magnetic field B at time t.
%  INPUT:
%   s1,s2: dmrg system
%   Phi:   state vectors
%   h:     n.n interaction between (function m-file)
%   M:     vector containing the dimensions of different blocks
%   B:  single-site time dep. magnetic field (function m-file).
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
  rhoL = PartialTrace(Phi*Phi',d*M(s1.l));
  [VL,FL] = eig(rhoL);
  [FL,p] = sort(-diag(FL));
  VL = VL(:,p);
  WL = VL(:,1:M(s1.l+1));
  WLM = kron(WL,speye(d));
  hLM = WLM'*kron(speye(M(s1.l)),feval(h,s1.l+1,s1.l+2))*WLM;
  if s1.l == s1.L-4
    hMR = feval(h,s1.l+3,s1.l+4);
  else
    WMR = kron(speye(d),s2.WR);
    hMR = WMR'*kron(feval(h,s1.l+3,s1.l+4),speye(M(s1.L-s1.l-4)))*WMR;
  end
  if nargin == 9
    hM = feval(h,s1.l+2,s1.l+3)+kron(feval(B,t,s1.l+2),speye(d))+...
         kron(speye(d),feval(B,t,s1.l+3));
    HL = WL'*(kron(s1.HL,speye(d))+s1.hLM+...
              kron(speye(M(s1.l)),feval(B,t,s1.l+1)))*WL;
  else
    hM = feval(h,s1.l+2,s1.l+3);
    HL = WL'*(kron(s1.HL,speye(d))+ s1.hLM)*WL;
  end
  r = dmrg_system(s1.L,s1.l+1,HL,hLM,hM,hMR,s2.HR,WL,s2.WR);


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

