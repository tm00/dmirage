function M = dmrg_dim(J,L,lm)
% DMRG_DIM - block dimensions for DMRG calculation
%   M = dmrg_dim(J,L,lm) sets the block dimensions for spin J, chain length
%   L, maximal full block length lm (lm does not have to be an integer)
  
  N = 2*J+1;
  
  % true max. full block length
  lt = floor(lm);
  % standard block dimension
  chi = floor(N^lm);
  M = repmat(chi, [L-2 1]);
  % smallest blocks have less total dimension than chi
  M(1:lt) = N.^(1:lt);
  % therefore larger blocks can have more states
  for l=lt:-1:1
    M(L-2-l) = min(floor(chi^2*N^(-l)),M(L-3-l)*N);% N.^(2*lm-(1:lm));
  end
  M(L-2) = M(L-3);
  
  % we want to keep total dimension constant
  for l=1:L-3
    M(l)*M(1)^2*M(L-l-2);
  end

