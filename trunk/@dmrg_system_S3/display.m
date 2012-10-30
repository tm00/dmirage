function display(s)
% DISPLAY - Command window display of DMRG system with conserved S3
%   
  
  display(s.dmrg_system)
  disp('left block total S3 =')
  disp(s.S3L)
  disp('right block total S3 =')
  disp(s.S3R)