function B = Bfield1(t,x)
  % XXZ parameters
  global J L bf v
  
  idL = speye(L);
  N = 2*J+1;
  a = 1:N-1;
  b = 2:N;
  Splus = sparse(a, b, sqrt((2*J-a+1).*a), N, N);
  Smin = (Splus).';
  S1=(1/2)*(Splus+Smin);
  
  if t <= v^(-1)
    B = idL(x,L/2)*(1-v*t)*bf*S1 + idL(x,L/2+1)*v*t*bf*S1;
  elseif t <= 2*v^(-1)
    B = idL(x,L/2+1)*(2-v*t)*bf*S1 + idL(x,L/2+2)*(v*t-1)*bf*S1;
  else
    error('Change file Bfield1.m to allow t>2*tau !')
  end
      
  