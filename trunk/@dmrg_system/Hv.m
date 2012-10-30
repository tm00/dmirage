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
