function [aa,vv] = nevis_initialize(b,s,gg,pp,oo)
% Assign default prescribed fields and boundary conditions and initial
% conditions [ probably need to change initial hs and phi ]
% Inputs
%   b bed elevation on nodes
%   s surface elevation on nodes
%   gg grid structure
%   pp parameters
%   oo options [optional]
% Outputs
%   aa structure containing prescribed fields and boundary conditions
%   vv structure containing initial variables
%
% 21 August 2014
if nargin<5, oo = struct; end
if ~isfield(pp,'u_b'), pp.u_b = 1; end
if ~isfield(pp,'l_c'), pp.l_c = 1; end
if ~isfield(pp,'phi_s'), pp.phi_s = -inf; end

%% reshape inputs
b = reshape(b,gg.nIJ,1);
s = reshape(s,gg.nIJ,1);

%% hydraulic potential
H = max(s-b,0);
phi_a = pp.r*b; 
phi_0 = phi_a+(s-b);

%% channel roughness sizes
lcx = pp.l_c*ones(gg.eIJ,1);
lcy = pp.l_c*ones(gg.fIJ,1);
lcs = pp.l_c*ones(gg.cIJ,1);
lcr = pp.l_c*ones(gg.cIJ,1);

%% storage capacity [ value taken from pp ]
sigma = pp.sigma*ones(gg.nIJ,1);
sigma(gg.nout) = NaN;

%% basal velocity [ value taken from pp ]
Ub = pp.u_b*ones(gg.nIJ,1); 
Ub(gg.nout) = 0;

%% basal melting rate [ value taken from pp ]
m = pp.melt*ones(gg.nIJ,1); 
m(gg.nout) = 0;

%% input
E = 0*ones(gg.nIJ,1); 
E(gg.nout) = 0;

%% variables
phi = phi_a+0*(phi_0-phi_a);
hs = 0*ones(gg.nIJ,1); 
Sx = 0*ones(gg.eIJ,1);
Sy = 0*ones(gg.fIJ,1);
Ss = 0*ones(gg.cIJ,1);
Sr = 0*ones(gg.cIJ,1);

%% boundary conditions
aa.phi = max(phi_a(gg.nbdy),pp.phi_s); 

%% put fields and variables in structs
aa.phi_0 = phi_0; 
aa.phi_a = phi_a; 
aa.lcx = lcx; 
aa.lcy = lcy; 
aa.lcs = lcs; 
aa.lcr = lcr; 
aa.b = b; 
aa.s = s; 
aa.H = H; 
aa.sigma = sigma; 
aa.m = m; 
aa.Ub = Ub; 
aa.E = E;

vv.phi = phi; 
vv.hs = hs; 
vv.Sx = Sx; 
vv.Sy = Sy; 
vv.Ss = Ss; 
vv.Sr = Sr;

end