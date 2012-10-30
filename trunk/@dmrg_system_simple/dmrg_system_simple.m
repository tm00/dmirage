function s = dmrg_system_simple(varargin)
% DMRG_SYSTEM_SIMPLE - Simple DMRG system constructor (no Hamiltonians)
%   dmrg_system(L,l,WL,WR) creates a dmrg system
%   consisting of:
%      L:   total chain length
%      l:   left block length
%      WL:  left block transformation matrix l -> l+1 ("WL{l}")
%      WR:  right block transformation matrix L-l-2 -> L-l-1
%           ("WR{L-2-l}") 
  
  switch nargin
   case 0 % if no input arguments, create empty system
    s.L = 0;
    s.l = 0;
    s.WL = [];
    s.WR = [];
    s = class(s,'dmrg_system_simple'); % create object
   case 1 % if single argument of class dmrg_system_simple, return it
    if (isa(varargin{1},'dmrg_system_simple'))
      s = varargin{1};
    else 
      error('Input argument is not a dmrg system')
    end
   case 4
    s.L = varargin{1};
    s.l = varargin{2};
    s.WL = varargin{3};
    s.WR = varargin{4};
    s = class(s,'dmrg_system_simple'); % create object
   otherwise
    error('Wrong number of input arguments')
  end
  
  
