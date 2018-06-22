%% Run file for Hewitt (2013) nevis model as used in Stevens et al., JGR, 2018.
% 140km-wide region with RACMO runoff input starting on DOY 1 2009.
% Starting condition is saved values from 2008 DOY 365 (2008 = a spin up year). 

format compact;
clear
oo.root = '';           % filename root
oo.fn = 'nevis_22222';  % filename
oo.code = './nevis';   % code directory
addpath(oo.code);

%% parameters
% default parameters
[pd,oo] = nevis_defaults([],oo);    
% alter default parmaeters 
pd.c_e_reg2 = 0.01/1e3/9.81;        % elastic sheet thickness [m/Pa]
pd.N_reg2 = 1e4; % 1e3              % regularisation pressure for elastic sheet thickness 
pd.u_b = 100/pd.ty;                 % sliding speed [m/s]
pd.sigma = 1e-3;                    % englacial void fraction
pd.h_r = 0.1;                       % roughness height [m]
pd.l_r = 10;                        % roughness length [m]
pd.l_c = 1000;                      % sheet width contributing to conduit melting [m] default = 10 m
pd.k_s = 1e-3;                      % sheet permeability constant
pd.tau_b = 60e3; %60 kPa            % driving stress [Pa]
pd.melt = pd.G/pd.rho_w/pd.L;       % geothermal heat derived basal melt [m/s]
pd.melt = (pd.G+(pd.u_b*pd.tau_b))/pd.rho_w/pd.L;  % geothermal heat + frictional heating derived basal melt [m/s]
pd.meltinterior = ((pd.G+((100/pd.ty)*pd.tau_b))/pd.rho_w/pd.L)*1e3; % flux of basal melt up to the ~icedivide (200 km) [m2/s]
% pd.gamma_cc = 7.8e-8;		        % beta, pressure dependent freezing pt.

% non-dimensionalise
[ps,pp] = nevis_nondimension(pd);   

%% grid and geometry
load('morlighem_for_nevis_140km'); % load Morlighem bedmap (previously collated)
dd = morlighem_for_nevis_140km; dd.skip = 6;
gg = nevis_grid(dd.X_km(1:dd.skip:end,1)/ps.x,dd.Y_km(1,1:dd.skip:end)/ps.x,oo); 
b = reshape(dd.B_km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);
s = reshape(dd.S_km(1:dd.skip:end,1:dd.skip:end)/ps.z,gg.nIJ,1);

%% mask with minimum ice thickness
H = max(s-b,0);
Hmin = 0/ps.z; 
nout = find(H<=Hmin);
gg = nevis_mask(gg,nout); 
gg.n1m = gg.n1;             % label all edge nodes as boundary nodes for pressure

%% label boundary nodes
gg = nevis_label(gg,gg.n1m);
oo.adjust_boundaries = 1;   % enable option of changing conditions

%% plot grid
%nevis_plot_grid(gg); return;  % check to see what grid looks like

%% initialize variables
[aa,vv] = nevis_initialize(b,s,gg,pp,oo);       % default initialisation
pd.k_f=0.9;                                     % percent overburden (k-factor) 
vv.phi = aa.phi_a+pd.k_f*(aa.phi_0-aa.phi_a);   % initial pressure  k_f*phi_0
N=aa.phi_0-vv.phi;                              % N for initial cavity sheet size 
vv.hs = ((((pd.u_b*pd.h_r/pd.l_r)./((pd.u_b/pd.l_r)+(pd.K_c.*((ps.phi*N).^3)))))./ps.h); % initial cavity sheet size as f(N)

%% boundary conditions
% aa.phi = aa.phi_a(gg.nbdy)+k_factor*(aa.phi_0(gg.nbdy)-aa.phi_a(gg.nbdy));    % prescribed boundary pressure
% aa.phi_b = aa.phi_0;             % prescribed boundary pressure at overburden
aa.phi_b = max(aa.phi_0,aa.phi_a); % prescribed boundary pressure at overburden or atmospheric (LAS 18 Nov. 2015)

%% moulins
oo.density_moulins = 1;
oo.keep_all_moulins = 0;
% create new moulin locations:
%[pp.ni_m,pp.sum_m] = nevis_moulins_density_joughin([],[],gg,oo,aa,ps);
% load consistent location moulins for all runs:
load('nevis/nevis_170207a.mat','pp');

%% surface input
% % to include distributed runoff from e.g. RACMO, define function
% % runoff(t,gg) to return runoff (m/s) at time t (s), at each point on the
% % grid (ie runoff(t,gg) should return a vector of size gg.nIJ-by-1), then
% % include as:
% pp.runoff_function = @(t) runoff(ps.t*t,gg)/ps.m; oo.runoff_function = 1;
load('nevis/runoff_2009_nevis140.mat');       % load data for year of interest (previously collated)
% RACMO distributed input
oo.distributed_input = 0;                     % If set to 1 turns on RACMO distributed input
pp.runoff_function = @(t) runoff(((t*ps.t)/pd.td),runoff_2009_nevis140)./ps.m;  % distributed input (m/sec)

% RACMO moulin input
oo.input_function = 1;                        % If set to 1 turns on RACMO moulin input (m3/sec)
pp.input_function = @(t) runoff_moulins(((t*ps.t)/pd.td),runoff_2009_nevis140,pp.sum_m,gg.Dx(1))./ps.m; % RACMO moulin input (m3/sec)

%% Timesteps and saving model output
oo.dt = 1/24*pd.td/ps.t; oo.save_timesteps = 1; oo.save_pts_all = 1; oo.pts_ni = pp.ni_m;

%% save initial parameters
save([oo.root,oo.fn],'pp','pd','ps','gg','aa','vv','oo');

%% timestep 
load('nevis/nevis_22221/0365.mat','vv','tt')
[tt,vv,info] = nevis_timesteps([1:1:365]*pd.td/ps.t,vv,aa,pp,gg,oo,pd,ps);

