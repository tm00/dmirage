function s = dmrg_system_S3(varargin)
% DMRG_SYSTEM_S3 - DMRG system with conserved S3 constructor
%   dmrg_system_S3(L,l,HL,hLM,hM,hMR,HR,WL,WR,S3L,S3R) creates a dmrg
%   system with conserved total S3 consisting of:
%      L:   total chain length
%      l:   left block length
%      HL:  left block Hamiltonian (size dl x dL)
%      hLM: interaction left block - middle sites (size d*dL x d*dL)
%      hM:  middle sites Hamiltonian (size d^2 x d^2)
%      hMR: interaction middle sites - right block (size d*dR x d*dR)
%      HR:  right block Hamiltonian (size dR x dR)
%      WL:  left block transformation matrix l -> l+1 ("WL{l}")
%      WR:  right block transformation matrix L-l-2 -> L-l-1
%           ("WR{L-2-l}") 
%      S3L: left block total S3
%      S3R: right block total S3
%      mag: magnetization sector
  
  switch nargin
   case 0 % if no input arguments, create empty system
    p = dmrg_system;
    s.S3L = [];
    s.S3R = [];
    s.mag = 0;
    s = class(s,'dmrg_system_S3',p);
   case 1 % if single argument of class dmrg_system_S3, return it
    if (isa(varargin{1},'dmrg_system_S3'))
      s = varargin{1};
    else 
      error('Input argument is not a dmrg system with conserved S3')
    end
   case 12
    p = dmrg_system(varargin{1:9});
    s.S3L = varargin{10};
    s.S3R = varargin{11};
    s.mag = varargin{12};
    s = class(s,'dmrg_system_S3',p);
   otherwise
    error('Wrong number of input arguments')
  end
