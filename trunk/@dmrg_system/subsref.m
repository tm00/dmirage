function val  = subsref(s, field)
% SUBSREF - Subscripted reference for dmrg system
%   
  switch field.type
   case '()'
    error('() not supported, use .')
   case '{}'
    error('{} not supported, use .')
   case '.'
    switch field.subs
     case 'L'
      val = s.L;
     case 'l'
      val = s.l;
     case 'HL'
      val = s.HL;
     case 'hLM'
      val = s.hLM;
     case 'hM'
      val = s.hM;
     case 'hMR'
      val = s.hMR;
     case 'HR'
      val = s.HR;
     case 'WL'
      val = s.WL;
     case 'WR'
      val = s.WR;
     otherwise
      error('Invalid field name')
    end
  end
  
  
