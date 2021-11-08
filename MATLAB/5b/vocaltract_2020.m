% Digital Waveguide Vocal Tract Model, based on the example by
% Mathur and Story in their paper:
%
% Mathur S., Story B. "Vocal tract modeling: Implementation of continuous
% length variations in a half-sample delay Kelly-Lochbaum model", ISSPIT 2003 
%
% Vocal tract section data from:
%
% Story, B., Titze, I.R., Parameterization of vocal tract area functions 
% by empirical orthogonal modes, Journal of Phonetics,26(3),1998, pp.
% 223-260, https://doi.org/10.1006/jpho.1998.0076.
%
% Adapted by DTM for LF excitation, Sound Output, Dynamics, Oct 2020
%

% Load crossectional area function, there are 44 sections in the profile.

AreaFile = load('areaFunctions.mat');
a = AreaFile.AreaFile.a;           % bart/father 
ae = AreaFile.AreaFile.ae;         % bat/lad
bird = AreaFile.AreaFile.bird;     % (3) bird
i_vowel = AreaFile.AreaFile.i;           % beet/see 
O = AreaFile.AreaFile.O;           % ball/law 
Q = AreaFile.AreaFile.Q;           % bod/not
u = AreaFile.AreaFile.u;           % food/soon
U = AreaFile.AreaFile.U;           % foot/put
V = AreaFile.AreaFile.V;           % but

%Define the waveguide parameters
Fs = 44100;                 % Sample Rate
rho = 1.2041;               % Density of Air at 20degC (kg/m^3)
c = 343.26;                 % Speed of sound at 20degC (m/s)
duration = 1.0;             % Sound output Duration (s)
nSegments = 44;             % Number of segments in the waveguide
nSamples = Fs*duration;     % No of speech samples (Pout) required

% Excitation
% The unit impulse/glottal volume velocity input
uin = zeros(nSamples, 1);
uin(1)= 200;
% Read in audio file for more advaned input
input_file = audioread('LFVib1000ms44100.wav'); % Use wave file as input (replaces uin)

%Upper and lower delay lines
u_delay = zeros(nSegments, 1); % Upper Right-going Delay Line
l_delay = zeros(nSegments, 1); % Lower Left-going Delay Line

Pout = zeros(nSamples, 1);     % Output

% Boundary Conditions
r_g = 0.99; % Reflection (with loss) at glottis end
r_l = -0.99; % Reflection (with loss) at lip end

% Define 2D matrix for cross-sectional area values
mesh_grid = zeros(nSamples, nSegments);
for n=1:nSegments
    mesh_grid(n) = (1-((n-1)/(nSegments-1)))*bird(n) + ((n-1)/(nSegments-1))*Q(n);
end    

%System Update Equations

for n=1:nSamples
    
    % Select a vowel for synthesis
    
    for j=1:nSegments
        A(j) = mesh_grid(j); 
    end
    
    % Reflection coefficients derived from the cross-sectional areas
    % The array of k-values are the reflection coefficients from tube section
    % to tube section as defined by the Kelly-Lochbaum scattering junction

    k = (A(1:(nSegments-1))-A(2:nSegments))./(A(1:nSegments-1)+A(2:nSegments));
    
    % Boundary Condition 1 - glottis end
    % Reflection from lower delay line to upper line at closed glottal end
    % + contribution from input u(n). 
    u_delay(1) = r_g*l_delay(1) + input_file(n)*(rho*c)/A(1); % Input file replaces uin here for more natural sound!

    % Scattering equations
    for i=1:nSegments-1
        temp  = k(i)*( u_delay(i)-l_delay(i+1) );
        u_delay(i+1) = u_delay(i) + temp;
        l_delay(i) = l_delay(i+1) + temp;
    end

    % Boundary Condition 2 - lip end
    % reflection from upper line to lower line at open lip
    l_delay(44) = r_l*u_delay(44);

    % Sum the upper and lower delay lines at the last section to get the output
    Pout(n) = u_delay(1) + l_delay(1);
    
end
Pout = Pout/(max(abs(Pout)));     % Normalise Pout

figure(1);
clf;
plot(Pout);
xlabel('time (samples)');
ylabel('Amplitude');

fftSize = 2^16;
f = (0:fftSize-1)*(Fs/fftSize);

figure(2);
clf;
plot(f,20*log10(abs(fft(Pout,fftSize))/max(abs(fft(Pout,fftSize)))));
axis([0 5000 -100 0]);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude Response (dB)');

soundsc(Pout, Fs);
audiowrite('2Dtest.wav', Pout, Fs);