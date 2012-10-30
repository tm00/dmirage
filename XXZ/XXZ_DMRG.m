function [s,Phi,E,t,Egs,S1prof,S2prof,S3prof,Eprof,norm_psi,overlap] = XXZ_DMRG(mag)
% XXZ DMRG Calculations

  global  J D L bf v

  D = 2; % anisotropy
  J = 0.5; % spin value
  L = 20; % chain length
  bf = 0.5; % strength of (moving) impurity

  sweep_steps = 3;
  time_steps = 200;
  ntrot = 10;
  v = 0.005; % velocity of moving impurity
  tf = 200;
  
  mag = 0.0;
  lm = 3;
  Ne = 3;
  
  M = dmrg_dim(J,L,lm);
  
  % display parameters
  disp(sprintf('\n%s\n','XXZ parameters:'))
  disp(['   L  = ' num2str(L)]);
  disp(['   J  = ' num2str(J)]);
  disp(['   D  = ' num2str(D)]);
  disp(['   bf = ' num2str(bf)]);
  disp(['   v  = ' num2str(v)]);
  disp(['   tf = ' num2str(tf)]);
  disp(sprintf('\n%s\n','DMRG parameters:'))
  disp(['   lm = '  num2str(lm)]);
  disp(['   Ne = '  num2str(Ne)]);
  disp(['   sweep_steps = '  num2str(sweep_steps)]);
  disp(['   time_steps  = '  num2str(time_steps)]);
  disp(['   ntrot       = '  num2str(ntrot)]);
  
  % first compute mag=0 ground state for system
  [s,Phi,E] = gsdmrgS3(M,Ne,'hXXZ',mag,sweep_steps,1);
  
  % time evolution
  [t,Egs,S1prof,S2prof,S3prof,Eprof,norm_psi,overlap] = ...
      tddmrg(M,Ne,s,'hXXZ','Bfield1',sweep_steps,tf,time_steps,ntrot);

  % save results
  save ( strcat('D=',num2str(D),'_L=',int2str(L),'_bf=',num2str(bf),'_v=', ...
                int2str(time_steps),'_lm=',int2str(lm),'_ntrot=', ...
                int2str(ntrot),'.mat') ) 
