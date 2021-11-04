function MassSpringDamper(F, z)
% A Single Mass-Spring-Damper Simple Harmonic Oscillator 
% DTM 1/10/2018

% System Definitions
% Fs = 44100;       % Sample Rate for audio output
Fs = F;
T = 1/Fs;           % Sample period
N = 44100;          % Number of samples in output
      
v1_init = 0;        % Initial Velocity of mass, v
x1_init = 1;        % Initial Displacement of the mass, x

% z1_init = 100;    % Initial damping of the mass, z
z1_init = z;
m1 = 30;            % Value of mass, m1

% Scaling parameters for mass and damping
z1_init = z1_init/1000000000;     % z given in N*s/(kg*10e-9)
m1 = m1/1000000000;               % mass given in kg*10e-9

% Calculate spring constant k
k = 0.5;
k = k*T*T;              % Scaling parameter according to FD expressions
z1_init = z1_init*T;    % Simplify the expression in the main loop

% Initialise displacement update variables for time t = 0 and t = -1
x1_1 = x1_init;        % Current displacement value, x1_1, for mass, m1
x1_2 = x1_init-(1000*v1_init*Fs); % Previous disp. value, x1_2, for m1 

% Initialise an array to capture output at a single point
out = zeros(1, N);

% Initialise an array for real-time graphing.
m1_out = zeros(1,250);

% Main Loop
for n=1:N 
    % Calculate total force acting on mass: F = F_{restoring} + F_{friction}
    % F = -kx + -zv
    % ENTER YOUR CODE HERE FOR CALCULATING f1, the force acting on m1
    v1_init = x1_1 - x1_2;
    f1 = (-k * x1_1) + (-z1_init * v1_init); 
    
    % Using the value for f1, calculate the next displacement for m1
    % ENTER YOUR CODE HERE FOR CALCULATING the mass displacement out(n)
    x1_next = (f1/m1) + 2*(x1_1) - x1_2;
    out(n) = x1_next;
    
    % Resamples the output at a rate of Fs/200 to plot the mass 
    % displacement in real-time.
    if mod(n,200)==0    
        m1_out(n/200) = x1_next; 
        plot(m1_out);
        drawnow;
    end
    
    % update values and prepare for next pass
    x1_2 = x1_1;
    x1_1 = x1_next;
    
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

