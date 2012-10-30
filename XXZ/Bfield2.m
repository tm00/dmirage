function B = Bfield2(t,x)
    
  global J L bf v
  
  N = 2*J+1;
  a = 1:N-1;
  b = 2:N;
  Splus = sparse(a, b, sqrt((2*J-a+1).*a), N, N);
  Smin = (Splus).';
  S1=(1/2)*(Splus+Smin);

  B=bf*0.25*(1+cos(pi*(x-v*t))).^2 .* chi1(x-v*t,L)*S1;
  
function y = chi1(x,L)
  y = zeros(size(x));
  y(L/2-1<=x & x<=L/2+1) = 1;
