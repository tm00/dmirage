function display(s)
% DISPLAY - Command window display of simple DMRG system
%   
  disp('SIMPLE DMRG SYSTEM:')
  disp('system length =')
  disp(s.L)
  disp('left block length =')
  disp(s.l)
  disp('left block transformation matrix =')
  disp(s.WL)
  disp('right block transformation matrix =')
  disp(s.WR)
