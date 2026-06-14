%% Physical constants
q    = 1.602176634e-19;      % C
kB   = 1.380649e-23;         % J/K
T    = 300;                  % K
VT   = kB*T/q;               % thermal voltage, V
vsat   = 1e5;                              % saturation velocity ([] = m/s)

eps0  = 8.8541878128e-12;    % F/m
epsSi = 11.7*eps0;           % F/m

gamma = 1;
tau_e = 0.2e-12; % [0.05 0.1 0.2 0.5 1.0]ps
TL = T;
%% Device geometry

small_device = 0; %set to 1 to simulate the small device

if small_device == 1
    L     = 1e-6;                % m
    LC    = 0.15e-6;             % m
    sigma = 10e-9;               % m
else 
    L     = 3e-6;
    LC    = 0.45e-6;
    sigma = 30e-9;
end

%% Doping
ND_low  = 1e22;              % m^-3 = 1e16 cm^-3
ND_high = 1e24;              % m^-3 = 1e18 cm^-3

%% Electron mobility and diffusion coefficient
mun_0 = 0.10;                  % m^2/(V s)
Dn  = mun_0*VT;                % m^2/s

%% Mesh
N = 1000;                     % total number of nodes
x = linspace(0,L,N).';
h = x(2)-x(1);
x_edge = 0.5*(x(1:end-1) + x(2:end));

Nint = N-2;                  % number of interior nodes

%% Smooth n+-n-n+ doping profile
SL = 0.5*(1 - tanh((x - LC)/sigma));
SR = 0.5*(1 + tanh((x - (L-LC))/sigma));

ND = ND_low + (ND_high - ND_low).*(SL + SR);

%% Boundary conditions for density
n_left  = ND_high;
n_right = ND_high;

%% Bias continuation values
Va_list = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5];  % V
%Va_list = linspace(0,0.5,20);
%Va_list = 0;         % use this for equilibrium only

%% Useful scales
LD_low  = sqrt(epsSi*VT/(q*ND_low));
LD_high = sqrt(epsSi*VT/(q*ND_high));
Vbi     = VT*log(ND_high/ND_low);
