% Toolbox D-MiRaGe
% Copyright (C) 2006 Tom Michoel, Bruno Nachtergaele and Wolfgang Spitzer
% Version 2006-10-19
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
%% This set of programs is free software; you can redistribute it and/or 
%% modify it under the terms of the GNU General Public License as
%% published by the Free Software Foundation; either version 2 of the
%% License, or (at your option) any later version.
%%
%% This set of programs is distributed in the hope that it will be
%% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License along
%% with this program; if not, write to the Free Software Foundation, Inc.,
%% 675 Mass Ave, Cambridge, MA 02139, USA.
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
%% You are welcome to use the code for your research under the terms of the
%% license. However, please acknowledge its use with the following
%% citation:
%%
%%   Tom Michoel, Bruno Nachtergaele and Wolfgang Spitzer; Transport of
%%      one-dimensional interfaces in the Heisenberg model;
%%      arXiv:ARXIV_NO (2007)
%%      (replace ARXIV_NO with the arXiv number where you downloaded this
%%       paper and software)
%%
%% Author contact information:
%%
%%  Tom Michoel
%%
%%  Bioinformatics & Evolutionary Genomics
%%  VIB/Ghent University
%%  Technologiepark 927
%%  B-9000 Gent, Belgium
%%
%%  Email: tom.michoel@psb.ugent.be
%%  URL:   http://www.psb.ugent.be/~tomic/
%%
%%
%%  Bruno Nachtergaele
%%
%%  Department of Mathematics
%%  University of California
%%  One Shields Avenue
%%  Davis, CA 95616-8633, USA
%%
%%  Email: bxn@math.ucdavis.edu
%%  URL:   http://math.ucdavis.edu/~bxn/
%%
%%
%%  Wolfgang Spitzer
%%
%%  Institut fur Theoretische Physik
%%  Universitat Erlangen-Nurnberg
%%  91058 Erlangen, Germany
%%
%%  Email: wolfgang.spitzer@physik.uni-erlangen.de
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% D-MiRaGe is a Matlab toolbox for performing both ground state and
% time-dependent Density Matrix Renormalization Group computations for
% one-dimensional quantum spin systems which need not be translation
% invariant. 
%
% 1. INSTALLATION
% ================
%
% To install, unpack in a directory DIR and add to Matlab search path or
% change into DIR and run: 
%
% >> Install
%
%
% 2. CLASSES
% ==========
%
% D-MiRaGe code is object-oriented. The basic class is a "dmrg_system"
% with fields:
%   L:   total chain length
%   l:   left block length
%   HL:  left block Hamiltonian
%   hLM: interaction left block - middle sites
%   hM:  middle sites Hamiltonian
%   hMR: interaction middle sites - right block
%   HR:  right block Hamiltonian
%   WL:  left block transformation matrix l -> l+1 ("WL{l}")
%   WR:  right block transformation matrix L-l-2 -> L-l-1
%           ("WR{L-2-l}")
% and functions:
%   compute_phi: compute target low-energy states of a dmrg system
%   enlarge:     insert 2 sites in the middle of a dmrg system
%   sweep_left:  perform one dmrg sweep step to the left
%   sweep_right: perform one dmrg sweep step to the right
%   Hv:          matrix-vector multiplication with the total Hamiltonian
%                without computing tensor products of block Hamiltonians
%
% For systems with symmetries, one should define derived classes of
% dmrg_system with overloading functions that exploit the particular
% symmetry. As an example, we provide "dmrg_system_S3" for a dmrg_system
% with conserved total S3-magnetization.
%
% Finally, there is a class "dmrg_system_simple" for dmrg_systems that no
% longer need to store block Hamiltonians, e.g., after a ground state
% computation to compute expectation values of one- and two-site operators.
%
% Class files are located in subdirectories whose name starts with '@'. For
% detailed usage information of any of the functions, type
%
% >> help function_name
%
% 
% 3. MAIN ALGORITHMS
% ==================
%
% The main commands for performing DMRG computations are:
%
% Ground state DMRG for n.n.-Hamiltonian without symmetries:
%
% >> [s,Phi,E] = gsdmrg(M,Ne,Hint,sweep_steps,pos);
%
% Ground state DMRG for n.n.-Hamiltonian which conserves total S3:
%
% >> [s,Phi,E] = gsdmrgS3(M,Ne,Hint,mag,sweep_steps,pos);
%
% Time-dependent DMRG for static n.n.-Hamiltonian with time-dependent
% external magnetic field, targeting the ground states of H(t) as
% explained in our paper:
%
% >> [s,Egs,S1prof,S2prof,S3prof,Eprof,norm_psi,overlap] = 
%        tddmrg(M,Ne,s,Hint,Bfield,sweep_steps,tf,time_steps,ntrot);
%
% Standard time-dependent DMRG for static n.n.-Hamiltonian with
% time-dependent external magnetic field, targeting the evolved state
% directly:
%
% >> [s,Sprof,Eprof] =
%        tddmrg_simple(M,Ne,s,Hint,Bfield,sweep_steps,tf,time_steps,ntrot);
%
% Ground state properties for a series of time-dependent Hamiltonians
% (static n.n. + time-dep field), using sweeping procedure as in tddmrg: 
%
% >> [s,energies,S3prof,enprof] =
%       gsdmrg_td(M,Ne,s,Hint,Bfield,sweep_steps,tf,time_steps,ntrot);
%
% There is one small helper function to compute block dimensions:
%
% >> M = dmrg_dim(J,L,lm);
%
% For details about input and output arguments of any of these functions,
% type:
%
% >> help function_name
%
%
% 4. MODEL SPECIFIC FUNCTIONS
% ===========================
%
% For actual model computations, a number of files need to be created,
% preferably stored in a separate directory not in the D-MiRaGe
% directory. As an example, we provide in the directory XXZ the files for
% the ferromagnetic, anisotropic Heisenberg Hamiltonian with kink boundary
% conditions and time-dependent magnetic fields as explained in our
% paper. Note that these functions use global model parameters.
%
% Compute nearest-neighbor-Hamiltonian (returns 0 matrix if y~=x+1):
%
% >> h = hXXZ(x,y);
%
% Compute single-site magnetic field perturbation at time t and position
% x (2 example fields provided): 
%
% >> B = Bfield1(t,x);
% >> B = Bfield2(t,x);
%
% Main program to execute DMRG computation for XXZ model:
%
% >> XXZ_DMRG();
%
% This program starts with a ground state computation to build an initial
% dmrg_system. This system is then passed as input to the tddmrg
% algorithm. The output is stored in a mat-file whose name is derived from
% the specific parameter settings.
%
% Finally, there is a small script that illustrates how to compute the
% total energy of the time evolved state from the energy and spin
% profiles outputted by tddmrg.
%

