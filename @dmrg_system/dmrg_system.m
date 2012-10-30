function s = dmrg_system(varargin)
% DMRG_SYSTEM - DMRG system constructor
%   dmrg_system(L,l,HL,hLM,hM,hMR,HR,WL,WR) creates a dmrg system
%   consisting of:
%      L:   total chain length
%      l:   left block length
%      HL:  left block Hamiltonian
%      hLM: interaction left block - middle sites
%      hM:  middle sites Hamiltonian
%      hMR: interaction middle sites - right block
%      HR:  right block Hamiltonian
%      WL:  left block transformation matrix l -> l+1 ("WL{l}")
%      WR:  right block transformation matrix L-l-2 -> L-l-1
%           ("WR{L-2-l}") 
  
  switch nargin
   case 0 % if no input arguments, create empty system
    s.L = 0;
    s.l = 0;
    s.HL = [];
    s.hLM = [];
    s.hM = [];
    s.hMR = [];
    s.HR = [];
    s.WL = [];
    s.WR = [];
    s = class(s,'dmrg_system'); % create object
   case 1 % if single argument of class dmrg_system, return it
    if (isa(varargin{1},'dmrg_system'))
      s = varargin{1};
    else 
      error('Input argument is not a dmrg system')
    end
   case 9
    p = dmrg_system_simple(varargin{1},varargin{2},varargin{8},varargin{9});
    s.L = varargin{1};
    s.l = varargin{2};
    s.HL = varargin{3};
    s.hLM = varargin{4};
    s.hM = varargin{5};
    s.hMR = varargin{6};
    s.HR = varargin{7};
    s.WL = varargin{8};
    s.WR = varargin{9};
    s = class(s,'dmrg_system'); % create object
   otherwise
    error('Wrong number of input arguments')
  end
  
  
