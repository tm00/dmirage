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
     case 'WL'
      val = s.WL;
     case 'WR'
      val = s.WR;
     otherwise
      error('Invalid field name')
    end
  end
  
  
