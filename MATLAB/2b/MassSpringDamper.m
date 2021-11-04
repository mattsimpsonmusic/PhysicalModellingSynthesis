function MassSpringDamper(F, z)
% A Single Mass-Spring-Damper Simple Harmonic Oscillator 
% Desirable inputs for testing are -> MassSpringDamper(44100, 100)

% System Definitions
Fs = F;             % Using input as sample rate
T = 1/Fs;           % Sample period
N = 44100;          % Number of samples in output
      
v1_init = 0;        % Initial Velocity of m1, v1
x1_init = 1;        % Initial Displacement of the m1, x1
v2_init = 0;        % Initial Velocity of m2, v2
x2_init = 1;        % Initial Displacement of m2, x2
    
z_init = z;         % Initial damping of the masses, z
m1 = 100;           % Value of first mass, m1
m2 = 50;           % Value of second mass, m2

% Scaling parameters for mass and damping
z_init = z_init/1000000000;       % z given in N*s/(kg*10e-9)
m1 = m1/1000000000;               % mass given in kg*10e-9
m2 = m2/1000000000;               % mass given in kg*10e-9

% Calculate spring constant k
k = 1;
k = k*T*T;            % Scaling parameter according to FD expressions
z_init = z_init*T;    % Simplify the expression in the main loop

% Initialise displacement update variables for time t = 0 and t = -1
x1_1 = x1_init;        % Current displacement value, x1_1, for mass, m1
x1_2 = x1_init-(1000*v1_init*Fs); % Previous disp. value for m1 
x2_1 = x2_init;
x2_2 = x2_init-(1000*v2_init*Fs); % Previous disp. value for m2 

% Initialise an array to capture output at a single point
out = zeros(1, N);

% Initialise arrays for real-time graphing.
m1_out = zeros(1,250);
m2_out = zeros(1,250);

% Main Loop
for n=1:N 
    % Calculate total force acting on masses: F = F_{restoring} + F_{friction}
    % F1 = -kx1 -zv1 + k(x2-x1) + z(v2-v1)
    % F2 = -k(x2-x1) - z(v2-v1)- kx2 - zv2
    % ENTER YOUR CODE HERE FOR CALCULATING f1, the force acting on m1
    v1_init = x1_1 - x1_2;
    v2_init = x2_1 - x2_2;
    
    f1 = (-k * x1_1) - (z_init * v1_init) + (k * (x2_1 - x1_1)) + (z_init * (v2_init - v1_init));
    f2 = (-k * (x2_1 - x1_1)) - (z_init * (v2_init - v1_init)) - (k * x2_1) - (z_init * v2_init);
    
    % Using the values for f1 & f2, calculate the next displacement for both m1 and m2
    % ENTER YOUR CODE HERE FOR CALCULATING the mass displacement out(n)
    x1_next = (f1/m1) + 2*(x1_1) - x1_2;
    x2_next = (f2/m2) + 2*(x2_1) - x2_2;
    
    out(n) = x2_next;
    
    % Resamples the output at a rate of Fs/200 to plot the mass 
    % displacement in real-time.
    if mod(n,200)==0    
        m1_out(n/200) = x1_next; 
        m2_out(n/200) = x2_next;
        plot(m1_out);
        hold on;
        plot(m2_out);
        drawnow;
        hold off;
    end
    
    % update values and prepare for next pass
    x1_2 = x1_1;
    x1_1 = x1_next;
    x2_2 = x2_1;
    x2_1 = x2_next;
    
end
% Plot time domain response of output
figure(2);
plot(out);
xlabel('Time (samples)');
ylabel('Amplitude');

fftSize = 8192;
f = (0:fftSize-1)*(Fs/fftSize);

% Plot frequency domain response of output
figure(3);
semilogx(f,20*log10(abs(fft(out,fftSize))/max(abs(fft(out,fftSize)))));
axis([0 8000 -40 0]);
xlabel('Frequency (Hz)');
ylabel('Magnitude Response (dB)');

sound(out, Fs);

end

