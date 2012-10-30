function [Phi, E] = compute_phi(varargin)
% COMPUTE_PHI - Compute target states of dmrg system with conserved S3
%   compute_phi(s,M,Ne,mag) : compute Ne target states and energies of system
%     s in magnetization sector mag
%   compute_phi(s,M,Ne,mag,opteigs) : compute Ne target states and energies
%     of system s in magnetization sector mag with opteigs options for eigs
%   compute_phi(s,M,Ne,mag,v0,opteigs) : compute Ne target states and energies
%     of system s in magnetization sector mag with trial vector v0 and
%     opteigs options for eigs  


  s = varargin{1};
  M = varargin{2};
  Ne = varargin{3};
  mag = varargin{4};
  switch nargin
   case 4
    opteigs = struct([]);
   case 5
    opteigs = varargin{5};
   case 6
    v0 = varargin{5};
    opteigs = varargin{6};
   otherwise
    error('Wrong number of input arguments');
  end
  
  t = s.dmrg_system;
  d = M(1);
  S3 = sparse(1:d,1:d,0.5*(d-1)-(0:d-1),d,d);
  S3M = kron(S3,eye(d)) + kron(eye(d),S3);

  if M(t.l) == 0 % all the way to the left there is no left block
    S3T =  kron(S3M,eye(M(t.L-t.l-2))) + kron(eye(d^2),s.S3R);
    H = sparse(kron(t.hM,eye(M(t.L-t.l-2))) + kron(eye(d),t.hMR) + ...
               kron(eye(d^2),t.HR));
  elseif M(t.L-t.l-2) == 0 % all the way to the right there is no right block
    S3T = kron(s.S3L,eye(d^2)) + kron(eye(M(t.l)),S3M);
    H = sparse(kron(t.HL,eye(d^2)) + kron(t.hLM,eye(d)) + ...
               kron(eye(M(t.l)),t.hM));
  else % they can't both be zero (L>=4)
    S3T = kron(s.S3L,eye(d^2*M(t.L-t.l-2))) + ...
          kron(eye(M(t.l)),kron(S3M,eye(M(t.L-t.l-2)))) ...
          + kron(eye(M(t.l)*d^2),s.S3R);
    H = sparse(kron(t.HL,eye(d^2*M(t.L-t.l-2))) + ...
               kron(t.hLM,eye(d*M(t.L-t.l-2))) + ...
               kron(eye(M(t.l)),kron(t.hM,eye(M(t.L-t.l-2)))) + ...
               kron(eye(M(t.l)*d),t.hMR) + kron(eye(M(t.l)*d^2), ...
						t.HR));
  end
  if nargin == 6
    [Phi,E] = jointEigs(H, S3T, mag, Ne, v0, opteigs);
  else
    [Phi,E] = jointEigs(H, S3T, mag, Ne, opteigs);
  end
  
function [Phi,E] = jointEigs(varargin)
% JOINTEIGS 
%   [Phi,E] = jointEigs(H, S3, mag, Ne, opteigs) compute for H the Ne
%       lowest eigenvalues in the sector 'mag' of S3
%   [Phi,E] = jointEigs(H, S3, mag, Ne, v0, opteigs) compute for H the Ne
%       lowest eigenvalues in the sector 'mag' of S3 using initial state
%       v0 
  
  H = varargin{1};
  S3 = varargin{2};
  mag = varargin{3};
  Ne = varargin{4};
  switch nargin
   case 5
    opteigs = varargin{5};
   case 6
    v0 = varargin{5};
    opteigs = varargin{6};
  end
  
  
  % S3 should be a diagonal matrix of half integers and H should commute
  % with S3 (both upto finite precision)
  [I,J] = find(S3);
  if I~=J
    error('S3 is not diagonal!')
  end
  % round values in S3 to the nearest half-integer, everything else is
  % finite precision errors
  s = full(diag(S3));
  s = round(s./0.5).*0.5;
  [I,J,K] = find(H>eps);
  % [H,S3] = 0 iff (s_i - s_j)H_(ij) = 0 for all i,j
  if (~isempty(find((s(I)-s(J)).*K > eps)))
    error('commutator of H and S3 is not zero!')
  end
  % we are in business now
  % find sector 'mag'
  ind = find(s == mag);
  % eigenvectors of S3 are given by standard basis
  Id = speye(size(S3));
  % initialize eigenvectors and eigenvalues of H
  Phi = zeros(length(H),Ne);
  E = zeros(Ne,1);
  % restrict H to sector and diagonalize
  Hs = sparse(Id(:,ind)'*H*Id(:,ind));
  % Hs could be non-symmetric due to finite precision
  Hs=0.5*(Hs+Hs');
  if nargin == 6
    opteigs.v0 = Id(:,ind)'*v0;
  end
  [V,D]=eigs(Hs,Ne,'sa',opteigs);
  Phi = Id(:,ind)*V;
  E = diag(D);
  
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
  % check the dimension of v
  if (length(v)~=dTot)
    error('input vector has wrong dimension');
  end
  % you'll need to play with tensor product indices to decipher this ...
  w = reshape(reshape(v, [dM*dR dL])*s.HL', [dL*dM*dR 1]) + ...
      reshape(reshape(v, [dMR dLM])*s.hLM', [dLM*dMR 1]) + ...
      reshape(reshape(s.hM*reshape(reshape(v, [dR dL*dM])', [dM dL*dR]), ...
                      [dL*dM dR])', [dL*dM*dR 1]) + ...
      reshape(s.hMR*reshape(v, [dMR dLM]), [dLM*dMR 1]) +  ...
      reshape(s.HR*reshape(v, [dR dL*dM]), [dL*dM*dR 1]);
