% linear code for benchmark


clear;  close all;


% Simulate mode selection
%  1: start from the very beginning of the simulation
%  2: continue from the last .mat file in ./  . Nothing should be 
%     changed except ndiagnose and visual
%  3: continue from the last .mat file in ../  . Used to simulate in 
%     different parameters
simulate_mode = 1;


% Automatic path management
if simulate_mode == 1
    if ~exist('code_linear/', 'file')
        error('No code_linear/ found. You might be in a further/ directory. Simulate_mode 1 unavailable.');
    end 
    code_path = fullfile(pwd, 'code_linear');
elseif simulate_mode == 2
    if ~exist('parameters.mat', 'file')
        error('No parameters.mat file found. Simulate_mode 2 unavailable.');
    end 
    load('parameters.mat', 'code_path', 'data_path');
elseif simulate_mode == 3
    if ~exist('../parameters.mat', 'file')
        error('No parameters.mat file found in ../  ,Simulate_mode 3 unavailable');
    end 
    load('../parameters.mat', 'code_path');
else
    error('Invalid simulate_mode');
end
addpath(code_path);

data_path = fullfile(pwd, 'data');
diagnose_path = fullfile(pwd, 'diagnose');
fig_path = fullfile(pwd, 'figs');
if ~exist(data_path, 'dir')
    mkdir('data')
end
if ~exist(diagnose_path, 'dir')
    mkdir('diagnose')
end
if ~exist(fig_path, 'dir')
    mkdir('figs')
end








% Mission
visual = 1;
enable_parallel = 1;


% Initial conditions
%  azimuthal mode number
m_number = 3;
%  parallel mode number
n_number = 1;
%  initial fluctuation mode amplitude
initamplitude_phi = 0.01;  % V
initamplitude_den = 0;


% Normalizing units
Tref = 3;  % reference temperature in eV
denref = 10e12;  % reference density in cm^-3
B0 = 1e3;  % Bz in Gauss
mu = 40;  % ion mass over proton mass. mu_argon=40


% Time step
dt = 0.004*2.611e-5;  % time step width in second
nt_per_diagnose = 300;
%  total diagnose times you want in current directory, which may not be
%  finished since the program might be interrupted
ndiagnose = 1800;


% Device settings
radius = 10;  % the radius of the vacuum vessel in cm
height = 270;  % the height of the cylinder vessel in cm


% Grid  (grids at x, y, and z boundaries are not included)
%  we use a cubic mesh, with x\in[-x_max, x_max]cm, y\in[-x_max, x_max]cm,
%  z\in[0,height]cm
x_max = 11;
%  number of grid points in both x and y directions(excludes boundaries)
%  has to be an even number
nx = 200; %146
nz = 22;


% Boundary conditions
%  x and y boundaries are given by the [calc] matrix
% z boundary conditions: 1 for free, 2 for periodic, 3 for sheath (not ready yet)
zbc_mode = 2;

% Diffusion
%  the Coulomb loarithm is always calculated using denref and Tref

%  dif_mode: model of perpendicular diffusion coefficient:
%  1: classical diffusion with den=0.5*denref and Te=0.5*Tref
%  2: updated using classical diffusion with local den and Te
%  3: two layers of constant diffusion coefficients:
%         dif_perp_in for r<rdif, dif_perp_out for r>rdif
dif_mode = 1;
rdif = 3.2;  % in cm
dif_perp_in = 6e3;  % in cm^2/s
dif_perp_out = 6e3;
% limit dif_perp
max_difperp = 6e3;  % in cm^2/s
min_difperp = 3e2;

%  parallel diffusion coefficient for den, w and vi
dif_z_in = 9.5e6;
dif_z_out = 9.5e6;%2e4;


% Conduction
%  parallel thermal conductivity \chi
rconduct = 3.2;  % in cm
conduct_z_in = 188.8888e6;  % in cm^2/s
conduct_z_out = 188.8888e6;


% Viscisity (uniform)
% viscosity is the perpendicular diffusion coefficient of vorticity and
%  ion parallel momentum given by ion-ion collisionality
viscosity = 4.8e3;  % in cm^2/s


% Decay (uniform)
%  loss of density by recombination
%  loss of temperature is unclear and is set to be at the same rate of density damping
%  loss of momentum by ion-neutral collsion

den_damp = 5e3;
momentum_damp = 3.4e3;


% Resistivity (by electron-ion collision \nu)
%  local_nustar=0: inv_nustar is calculated using denref and Tref
%  local_nustar=1: inv_nustar is calculated using local den and Te
%   Coulomb logarithm is always calculated using denref and Tref
local_nustar = 1;

main;
