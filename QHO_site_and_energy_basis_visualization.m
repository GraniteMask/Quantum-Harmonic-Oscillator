%% Constants
hbar = 6.58E-16; % [eV*s] the reduced Planck constant
c_cm = 2.997E10; % [cm/s] the speed of light
c_nm = 2.997E17; % [nm/s] the speed of light
amu2eV = 931.5E6; % [eV] convert amu to eV [E = m * c^2] 

%% Parameters
N = 71; % [] number of basis states in a truncated representation of QHO
m_amu = 6; % [amu] oscillator mass in amu
m = m_amu * amu2eV / (c_nm^2); % [eV*s^2/nm^2] oscillator mass

f_inv_cm = 300; % [cm^{-1}] oscillator frequency (f = c/lambda)
f = f_inv_cm * c_cm; % [Hz] oscillator frequency
w = 2*pi*f;  % [rad/s] oscillator angular frequency

nt = 75; % [] number of time points
TimeSpan = 5; % [To] calculation time span in units of To

%% Calculations

% Construct operators for the QHO
a = diag(sqrt(1:N-1), 1); % [] lowering operator in energy basis
ad = a';  % [] raising operator in energy basis

H = hbar * w * (ad * a + 0.5*eye(N)); % [eV] Hamiltonian (energy basis)

X = sqrt(hbar/(2*m*w))*(ad + a); % [nm] position operator (energy basis)

alpha = 1 + 1i*0;

% UEQ - eigenvectors of the X operator (energy basis)
[UEX, Xeigval] = eig(X);
UXE = UEX'; % transformation matrix from energy to site basis
% To transform psi_E in energy basis to psi_X in position basis, use
% psi_E = UXE * psi_X

Eo = hbar*w * 0.5; % [eV] energy scale for the problem
To = hbar/Eo; % [s] time scale for the problem


n = 2; % index for an eigenstate
eigenstate_n_E = zeros(N, 1);
eigenstate_n_E(n+1) = 1;
pd_eigenstate_n_E = conj(eigenstate_n_E) .* eigenstate_n_E;

% transform state from energy basis to position space
eigenstate_n_X = UXE * eigenstate_n_E;
pd_eigenstate_n_X = conj(eigenstate_n_X) .* eigenstate_n_X;


t = linspace(0, TimeSpan*To, nt); % [s] time vector

% Preparing the initial state
psi_0_case = 4;
psi_0 = zeros(N, 1); % zero vector for initial state
switch psi_0_case
    case 1 % initial state is the ground state
        psi_0(1) = 1;
    case 2 % initial state is the first excited state
        psi_0(2) = 1; 
    case 3 % linear combination of stationary state
        state_A = psi_0;
        state_B = psi_0;
        state_A(1) = 1;
        state_A(2) = 1;
        psi_0 = 1/sqrt(2) * (state_A + state_B);
    case 4 % Glauber state
        psi_gnd = psi_0; % ground state
        psi_gnd(1) = 1;
        % Displacement operator
        D = expm(alpha * ad - conj(alpha)*a);
        psi_0 = D * psi_gnd;
        
end

psi_t_E = zeros(N, nt); % time-varying state psi(t) [energy basis]
psi_t_X = zeros(N, nt); % time-varying state psi(t) [site basis]
for t_idx = 1:nt
    Ut = expm(-1i * H * t(t_idx)/hbar );
    % calculate t_idx time step [energy basis]
    psi_t_E(:, t_idx) = Ut * psi_0;
    % calculate t_idx time step [site basis]
    psi_t_X(:, t_idx) = UXE * psi_t_E(:, t_idx);
end

% probability densities
pd_t_E = conj(psi_t_E) .* psi_t_E;
pd_t_X = conj(psi_t_X) .* psi_t_X;


%% Visualization

% Energy basis visualization
subplot(1,2,1);
%bar3color(pd_t_E);
bar3(pd_t_E);
grid on;
set(gca, 'FontSize', 16, 'FontName', 'Times');
xlabel('Time step', 'Interpreter', 'latex');
ylabel('$n$', 'Interpreter', 'latex');

% Site basis visualization
subplot(1,2, 2);
pcolor(diag(Xeigval)*ones(1, nt), ...
    ones(N,1)* (t/To), ...
    pd_t_X);
set(gca, 'FontSize', 16, 'FontName', 'Times');
shading interp;
% grid on;
xlim(0.05*[-1,1])
ylabel('$t/T_o$', 'Interpreter', 'latex');
xlabel('$x$ (nm)', 'Interpreter', 'latex');
zlabel('Probability', 'Interpreter', 'latex');