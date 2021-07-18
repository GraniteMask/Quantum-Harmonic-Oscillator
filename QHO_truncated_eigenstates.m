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

%% Calculations

% Construct operators for the QHO
a = diag(sqrt(1:N-1), 1); % [] lowering operator in energy basis
ad = a';  % [] raising operator in energy basis

H = hbar * w * (ad * a + 0.5*eye(N)); % [eV] Hamiltonian (energy basis)

X = sqrt(hbar/(2*m*w))*(ad + a); % [nm] position operator (energy basis)

% UEQ - eigenvectors of the X operator (energy basis)
[UEX, Xeigval] = eig(X);
UXE = UEX'; % transformation matrix from energy to site basis
% To transform psi_E in energy basis to psi_X in position basis, use
% psi_E = UXE * psi_X


n = 2; % index for an eigenstate
eigenstate_n_E = zeros(N, 1);
eigenstate_n_E(n+1) = 1;
pd_eigenstate_n_E = conj(eigenstate_n_E) .* eigenstate_n_E;

% transform state from energy basis to position space
eigenstate_n_X = UXE * eigenstate_n_E;
pd_eigenstate_n_X = conj(eigenstate_n_X) .* eigenstate_n_X;


%% Visualization

% Energy basis visualization
subplot(1,3,1);
barh(0:N-1, pd_eigenstate_n_E);
grid on;
set(gca, 'FontSize', 16, 'FontName', 'Times');
ylabel('$n$', 'Interpreter', 'latex');
xlabel('Probability', 'Interpreter', 'latex');

% Site basis visualization
subplot(1,3, 2:3);
bar(diag(Xeigval), pd_eigenstate_n_X);
set(gca, 'FontSize', 16, 'FontName', 'Times');
grid on;
xlabel('$X$ (nm)', 'Interpreter', 'latex');
ylabel('Probability', 'Interpreter', 'latex');