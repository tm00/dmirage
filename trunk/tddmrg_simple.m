function [s,Sprof,Eprof] = tddmrg_simple(varargin)
% TDDMRG_SIMPLE - Simple Time Dependent Density Matrix Renormalization Group
%   TDDMRG(M,Ne,s,Hint,Bfield,sweep_steps,tf,time_steps,ntrot)
%     M: vector with dimensions for each possible block length;
%        note that spin 2*J+1=M(1) and system length L =
%        length(M)+2.
%     Ne: number of targeted states.
%     s: cell array of dmrg structures, such as computed by GSDMRG.
%     Hint: nearest neighbor interaction H_{x,x+1} (function
%           handle)
%     Bfield: time dependent external magnetic field
%     sweep_steps: number of sweeps through the system
%     tf: final time
%     time_steps: number of time steps between 0 and tf
%     ntrot: number of Trotter steps in each time step
% 
%   TDDMRG(M,Ne,s,Hint,Bfield1,Bfield2,sweep_steps,tf,time_steps,ntrot)
%     allows for the possibility to use different magnetic fields in the
%     sweeping and Trotter loops; Bfield1 for sweeping and Bfield2 for
%     Trotter 
%
%   OUTPUT:
%     we can output many things; for now we return the updated cell array
%     of dmrg structures and the S3 profile. 
  
  disp(sprintf('\n%s\n','Start time-dependent DMRG computation'))
  
  M = varargin{1};
  Ne = varargin{2};
  s = varargin{3};
  Hint = varargin{4};
  Bfield1 = varargin{5};
  switch nargin
   case 9
    Bfield2 = varargin{5};
    sweep_steps = varargin{6};
    tf = varargin{7};
    time_steps = varargin{8};
    ntrot = varargin{9};
   case 10
    Bfield2 = varargin{6};
    sweep_steps = varargin{7};
    tf = varargin{8};
    time_steps = varargin{9};
    ntrot = varargin{10};  
   otherwise
    error('Wrong number of input arguments')
  end
  N = M(1); % dimension of single-site Hilbert space
  J = 0.5*(N-1);
  a = 1:N-1;
  b = 2:N;
  Splus = sparse(a, b, sqrt((2*J-a+1).*a), N, N);
  Smin = (Splus).';
  S1 = 0.5*(Splus+Smin);
  S2 = -0.5*i*(Splus-Smin);
  S3 = sparse(1:N,1:N,J-(0:2*J),N,N);
  L = length(M)+2; % total system size;
  delta = tf/time_steps; % length of a time step
  dtrot = delta/ntrot; % length of trotter step;
  opteigs.disp = 0; % turn off eigs output
  opteigs.sigma = 'sr'; % which eigenvalues to compute
% $$$   S3prof = zeros(L,time_steps+1);   % store the S3 profile
% $$$   norm_psi = ones(1,time_steps+1); % store the norm of the evolved state
% $$$   overlap = ones(1,time_steps+1); % store the overlap with the ground state
  
  % the input dmrg system could be a ground state of H and not of
  % H(0)=H+B(0); therefore we start with sweeps to get the initial state;
  disp('-> sweep to initial state')
  if (isa(s{1},'dmrg_system_S3'))
    [Phi,E] = compute_phi(s{1},M,Ne,s{1}.mag,opteigs);
    s{1} = s{1}.dmrg_system;
  else
    [Phi,E] = compute_phi(s{1},M,Ne,opteigs);
  end
  for st = 1:sweep_steps
    for l = 1:L-4
      s{l+1} = sweep_right(s{l},s{l+1},Phi,Hint,M,Ne,opteigs,Bfield1,0);
      opteigs.v0 = AL1AR(s{l+1}.WL',s{l}.WR,N,Phi(:,1));
      [Phi,E] = compute_phi(s{l+1},M,Ne,opteigs);
    end
    for l = L-3:-1:2
      s{l-1} = sweep_left(s{l},s{l-1},Phi,Hint,M,Ne,opteigs,Bfield1,0);
      opteigs.v0 = AL1AR(s{l}.WL,s{l-1}.WR',N,Phi(:,1));
      [Phi,E] = compute_phi(s{l-1},M,Ne,opteigs);
    end
  end
  
  % in the remainder we need no more block Hamiltonians -> convert to
  % simple dmrg system
  for l=1:length(s)
    s{l} = dmrg_system_simple(s{l}.L, s{l}.l, s{l}.WL, s{l}.WR);
  end
  % initial condition
  Psi = Phi;
  % its S3 profile
%  S3prof(:,1) = exp_profile(S3,Psi(:,1),s,M);
  disp('-> start time evolution')
  for n=1:time_steps
    disp(['     (time step ', int2str(n), ')'])
    % time
    t = n*delta;
    
    % we update the block bases to the current time evolved states
    for st=1:sweep_steps
      for l=1:L-4
        s{l+1} = sweep_right(s{l},s{l+1},Psi,M);
        % update Psi
	Psi = AL1AR(s{l+1}.WL',s{l}.WR,N,Psi);
        for k=1:size(Psi,2)
	  %Psi(:,k) = AL1AR(s{l+1}.WL',s{l}.WR,N,Psi(:,k));
	  Psi(:,k) = Psi(:,k)/norm(Psi(:,k));
	end
      end
      for l=L-3:-1:2
        s{l-1} = sweep_left(s{l},s{l-1},Psi,M);
	% update Psi
	Psi = AL1AR(s{l}.WL,s{l-1}.WR',N,Psi);
        for k=1:size(Psi,2)
	  %Psi(:,k) = AL1AR(s{l}.WL,s{l-1}.WR',N,Psi(:,k));
	  Psi(:,k) = Psi(:,k)/norm(Psi(:,k));
	end
      end
    end
    
    % apply e^(-i*delta*H(t)) to Psi(t-delta) to get Psi(t); this is done by
    % dividing delta into ntrot smaller steps; notice this uses Bfield2
    for nt = 1:ntrot
      % compute exp(-i*dtrot*H)*Psi, ntrot times
      % apply propagator by sweeping throught the lattice
      % first links (1,2) and (2,3) can be applied in the same basis
      hM = feval(Hint,1,2) + kron(feval(Bfield2,t,1),eye(N));
      Psi = A1v(expm(-0.5*i*dtrot*hM),Psi,M);
      hM = feval(Hint,2,3) + kron(feval(Bfield2,t,2),eye(N));
      Psi = Alv(expm(-0.5*i*dtrot*hM),Psi,1,M);
      % sweep right
      for l=2:L-3
        % transform Psi
	Psi = AL1AR(s{l}.WL',s{l-1}.WR,N,Psi);
        hM = feval(Hint,l+1,l+2) + kron(feval(Bfield2,t,l+1),eye(N));
        % propagate Psi
	Psi = Alv(expm(-0.5*i*dtrot*hM),Psi,l,M);
      end
      % finally link (L-1,L), in the same basis as the last step
      hM = feval(Hint,L-1,L) + kron(feval(Bfield2,t,L-1),eye(N)) + ...
           kron(eye(N),feval(Bfield2,t,L));
      Psi = ALv(expm(-i*dtrot*hM),Psi,M);
      % (L-2,L-1) still in the same basis
      hM = feval(Hint,L-2,L-1) + kron(feval(Bfield2,t,L-2),eye(N));
      Psi = Alv(expm(-0.5*i*dtrot*hM),Psi,L-3,M);
      % sweep left
      for l=L-4:-1:1
        % basis transformation
	Psi = AL1AR(s{l+1}.WL,s{l}.WR',N,Psi);
	%Psi = kron(s{l+1}.WL,kron(eye(N),s{l}.WR'))*Psi;
	hM = feval(Hint,l+1,l+2) + kron(feval(Bfield2,t,l+1),eye(N));
        % propagate Psi
	Psi(:,k) = Alv(expm(-0.5*i*dtrot*hM),Psi(:,k),l,M);
	%Psi = kron(speye(M(l)), kron(expm(-0.5*i*dtrot*hM), speye(M(L-l-2))))*Psi;
      end
      % finally link (1,2)
      hM = feval(Hint,1,2) + kron(feval(Bfield2,t,1),eye(N));
      Psi = A1v(expm(-0.5*i*dtrot*hM),Psi,M);
      %Psi = kron(expm(-0.5*i*dtrot*hM), speye(N*M(L-3)))*Psi;
    end
    % compute its profile
    %S3prof(:,n+1) = exp_profile(S3,Psi(:,1),s,M);
  end
  
  % compute profile, energy of final state
  S1prof = exp_profile(S1,Psi(:,1),s,M);
  S2prof = exp_profile(S2,Psi(:,1),s,M);
  S3prof = exp_profile(S3,Psi(:,1),s,M);
  Sprof = [S1prof, S2prof, S3prof];
  h = feval(Hint,1,2); % Hint is translation invariant
  Eprof = exp2_profile(h,Psi(:,1),s,M);

  disp(sprintf('\n%s\n', 'End'))

function prof = exp_profile(Op,Psi,s,M)
% EXP_PROFILE - Compute expectation profile
%   EXP_PROFILE(Op,Psi,s,M) computes the expectation profile of the single
%   site operator Op in the state Psi, using the structure of dmrg_system
%   s with block dimensions M.

  L = length(M)+2; % system length
  d = M(1); % single site dimension
  prof = zeros(L,1);
  Psi1 = Psi;
  
  prof(1) = real((Psi1'*kron(Op,eye(d^2*M(L-3)))*Psi1)/norm(Psi1)^2);
  prof(2) = real((Psi1'*kron(eye(d),kron(Op,eye(d*M(L-3))))*Psi1)/ ...
                 norm(Psi1)^2);
  for l=2:L-3
    % basis transformation
    Psi1 = AL1AR(s{l}.WL',s{l-1}.WR,d,Psi1);
    %Psi1 = kron(s{l}.WL',kron(eye(d),s{l-1}.WR))*Psi1;
    % Op(l+1)
    prof(l+1) = real((Psi1'*kron(eye(M(l)),kron(Op,eye(d*M(L-2-l))))* ...
                      Psi1)/norm(Psi1)^2);
  end
  prof(L-1) = real((Psi1'*kron(eye(M(L-3)*d),kron(Op,eye(d)))*Psi1)/ ...
                   norm(Psi1)^2);
  prof(L) = real((Psi1'*kron(eye(M(L-3)*d^2),Op)*Psi1)/ ...
                 norm(Psi1)^2);

function prof = exp2_profile(Op,Psi,s,M)
% EXP2_PROFILE - Compute expectation profile
%   EXP2_PROFILE(Op,Psi,s,M) computes the expectation profile of the two
%   site operator Op in the state Psi, using the structure of dmrg_system
%   s with block dimensions M.

  L = length(M)+2; % system length
  d = M(1); % single site dimension
  prof = zeros(L-1,1);
  Psi1 = Psi;
  
  prof(1) = real(Psi1'*A1v(Op,Psi1,M))/norm(Psi1)^2;
  prof(2) = real(Psi1'*Alv(Op,Psi1,1,M))/norm(Psi1)^2;
  for l=2:L-3
    % basis transformation
    Psi1 = AL1AR(s{l}.WL',s{l-1}.WR,d,Psi1);
    prof(l+1) = real(Psi1'*Alv(Op,Psi1,l,M))/norm(Psi1)^2;
  end
  prof(L-1) = real(Psi1'*ALv(Op,Psi1,M))/norm(Psi1)^2;

  
function w = AL1AR (AL, AR, N, v)
% compute w = kron(AL,kron(eye(N),AR))*v
% by computing w = kron(AL,eye())*kron(eye(),AR)*v
  
  dL1 = size(AL,1);
  dL2 = size(AL,2);
  dR1 = size(AR,1);
  dR2 = size(AR,2);
  Ne = size(v,2); % number of states

  w = zeros(dL1*N*dR1, Ne);
  for k=1:Ne
    w0 = AR*reshape(v(:,k), [dR2 dL2*N]);
    w1 = reshape(w0, [dL2*N*dR1 1]);
    w2 = reshape(w1, [N*dR1 dL2])*AL.';
    w(:,k) = reshape(w2, [dL1*N*dR1 1]);
  end
  
  
function w = A1v (A,v,M)
% matrix vector multiplication with A_{1,2} \otimes 1

  L = length(M)+2;
  dM = M(1)^2;
  dR = M(1)*M(L-3);
  Ne = size(v,2); % number of states
  
  w = zeros(dM*dR,Ne);
  for k=1:Ne
    w0 = reshape(v(:,k), [dR dM])*A.';
    w(:,k) = reshape(w0, [dM*dR 1]);
  end

  
function w = Alv (A,v,l,M)
% matrix vector multiplication with 1 \otimes A_{l+1,l+2} \otimes 1

  L = length(M)+2;
  dM = M(1)^2;
  dL = M(l);
  dR = M(L-l-2);
  Ne = size(v,2); % number of states
  
  w = zeros(dL*dM*dR,Ne);
  for k=1:Ne
    w0 = reshape(v(:,k), [dR dL*dM]).';
    w1 = A*reshape(w0, [dM dL*dR]);
    w2 = reshape(w1, [dL*dM dR]).';
    w(:,k) = reshape(w2, [dL*dM*dR 1]);
  end
  
function w = ALv (A,v,M)
% matrix vector multiplication with 1 \otimes A_{L-1,L} 

  L = length(M)+2;
  dM = M(1)^2;
  dL = M(L-3)*M(1);
  Ne = size(v,2); % number of states
  
  w = zeros(dL*dM,Ne);
  for k=1:Ne
    w0 = A*reshape(v(:,k), [dM dL]);
    w(:,k) = reshape(w0, [dL*dM 1]);
  end
