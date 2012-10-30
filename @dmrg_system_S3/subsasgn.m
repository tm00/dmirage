function s  = subsasgn(s, field, value)
% SUBSASGN - Subscripted asignment for magnetization of dmrg system with conserved S3
%   

  switch field.type
   case '()'
    error('() not supported, use .')
   case '{}'
    error('{} not supported, use .')
   case '.'
    switch field.subs
     case 'mag'
      s.mag = value;
     otherwise
      error('Invalid field name')
    end
  end
  
  
