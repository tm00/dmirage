function val  = subsref(s, field)
% SUBSREF - Subscripted reference for dmrg system with conserved S3
%   
  u = s.dmrg_system;
  switch field.type
   case '()'
    error('() not supported, use .')
   case '{}'
    error('{} not supported, use .')
   case '.'
    switch field.subs
     case 'dmrg_system'
      val = s.dmrg_system;
     case 'L'
      val = u.L;
     case 'l'
      val = u.l;
     case 'HL'
      val = u.HL;
     case 'hLM'
      val = u.hLM;
     case 'hM'
      val = u.hM;
     case 'hMR'
      val = u.hMR;
     case 'HR'
      val = u.HR;
     case 'WL'
      val = u.WL;
     case 'WR'
      val = u.WR;
     case 'S3L'
      val = s.S3L;
     case 'S3R'
      val = s.S3R;
     case 'mag'
      val = s.mag;
     otherwise
      error('Invalid field name')
    end
  end
  
  
