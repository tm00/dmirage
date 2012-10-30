function [s,Phi,E] = gsdmrgS3(varargin)
% GSDMRGS3 - Ground State Density Matrix Renormalization Group with fixed S3
%   [s,Phi,E] = GSDMRGS3(M,Ne,Hint,mag,sweep_steps,pos) 
%      Compute Ne DMRG states Phi with energies E for the Hamiltonian
%      with n.n.-interaction Hint
%      INPUT:
%        M: vector with dimensions for each possible block length;
%           note that spin 2*J+1=M(1) and system length L =
%           length(M)+2.
%        Ne: number of targeted states.
%        Hint: nearest neighbor interaction H_{x,x+1} (function
%              handle)
%        mag: magnetization sector of S3 total to consider
%        sweep_steps: number of sweeps through the system
%        pos: end position of left block in which final configuration
%             should be (1,...,L-3).
%
%      OUTPUT: 
%        s: a cell-array of dmrg structures where s{l} corresponds
%           to a system with left block length l.
%        Phi, E: targeted states and energies in final configuration


  disp(sprintf('\n%s\n','Start ground state DMRG computation'))
  
  M = varargin{1};
  Ne = varargin{2};
  Hint = varargin{3};
  mag = varargin{4};
  sweep_steps = varargin{5};
  pos = varargin{6};  

  N = M(1); % dimension of single-site Hilbert space
  J = 0.5*(N-1);
  S3 = sparse(1:N,1:N,J-(0:2*J),N,N);
  L = length(M)+2; % total system size; has to be even and greater than 8
  opteigs.disp = 0; % turn off eigs output
  opteigs.sigma = 'SA'; % which eigenvalues to compute
 
  
  % no matter what 'mag' is asked, we start by building the ground state in sector
  % mag=0 (because not all allowed mag in system of length L are available in smaller systems)
  s{1} = dmrg_system_S3(4,1,zeros(M(1)),feval(Hint,1,2),feval(Hint,2,3),...
                        feval(Hint,3,4),zeros(M(1)),eye(M(1)),eye(M(1)),...
                        S3,S3,0.0);
 
  % enlarge system from 4 sites to L sites
  disp(['-> enlarging from 4 sites to ', int2str(L), ' sites'])
  disp('     (mag = 0)')
  for l=1:L/2-2
    s{l+1} = enlarge(s{l},Hint,M,Ne,opteigs);
    s{l+1}.mag = s{l}.mag;
  end
  [Phi,E] = compute_phi(s{L/2-1},M,Ne,s{L/2-1}.mag,opteigs);
  % now we can switch to desired mag value
  s{L/2-1}.mag = mag;
  % initialization sweep
  disp('-> initialization sweep')
  if (mag == round(mag))
    disp(['     (mag = ', int2str(mag), ')'])
  else
    disp(['     (mag = ', int2str(2*mag), '/2)'])
  end
  % don't use initial vectors v0 right after switching magnetizations
  for l = L/2-1:L-4
    s{l+1} = sweep_right(s{l},s{L-l-3},Phi,Hint,M,Ne,opteigs);
    [Phi,E] = compute_phi(s{l+1},M,Ne,s{l+1}.mag,opteigs);
  end
  for l = L-3:-1:2
    s{l-1} = sweep_left(s{l},s{l-1},Phi,Hint,M,Ne,opteigs);
    [Phi,E] = compute_phi(s{l-1},M,Ne,s{l-1}.mag,opteigs);
  end
  % further sweeps
  disp('-> convergence sweeps')
  % by now it should be safe to start using initial vectors v0
  for st = 1:sweep_steps
    disp(['     (sweep ', int2str(st), ')'])
    for l = 1:L-4
      s{l+1} = sweep_right(s{l},s{l+1},Phi,Hint,M,Ne,opteigs);
      v0 = kron(s{l+1}.WL',kron(eye(N),s{l}.WR))*Phi(:,1);
      [Phi,E] = compute_phi(s{l+1},M,Ne,s{l+1}.mag,v0,opteigs);
    end
    for l = L-3:-1:2
      s{l-1} = sweep_left(s{l},s{l-1},Phi,Hint,M,Ne,opteigs);
      v0 = kron(s{l}.WL,kron(eye(N),s{l-1}.WR'))*Phi(:,1);
      [Phi,E] = compute_phi(s{l-1},M,Ne,s{l-1}.mag,v0,opteigs);
    end
  end
  % sweep to desired end position
  disp(['     (sweep to position ', int2str(pos), ')'])
  for l = 1:pos-1
    s{l+1} = sweep_right(s{l},s{l+1},Phi,Hint,M,Ne,opteigs);
    v0 = kron(s{l+1}.WL',kron(eye(N),s{l}.WR))*Phi(:,1);
    [Phi,E] = compute_phi(s{l+1},M,Ne,s{l+1}.mag,v0,opteigs);
  end

  disp(sprintf('\n%s\n', 'End'))
