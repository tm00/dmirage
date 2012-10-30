function h = hXXZ(x,y)
% HXXZ - XXZ n.n. interaction
%   hXXZ(x,y=x+1) = h_{x,x+1}
  
  % XXZ parameters
  global J D
  N = 2*J+1; % dimension of single-site Hilbert space
  % Single-site spin matrices:
  S3 = sparse(1:N,1:N,J-(0:2*J),N,N);
  a = 1:N-1;
  b = 2:N;
  Splus = sparse(a, b, sqrt((2*J-a+1).*a), N, N);
  Smin = (Splus).';
  
  % twisted boundary condition
  bc =  J*sqrt(1-1/D^2)*(kron(S3,eye(N)) - kron(eye(N),S3)); 

  % nearest neighbor interaction
  if y == x+1
    h = -(1/(2*D))*( kron(Splus,Smin) + kron(Smin,Splus) ) - ...
        kron(S3,S3) + J^2*eye(N^2) + bc;
  else 
    h == zeros(N^2);
  end
