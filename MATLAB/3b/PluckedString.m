% Plucked String Digital Waveguide Sound Synthesis
% Based on travelling wave, two delay-line model
% Variable plucking position and output possible
% Simple bridge reflection coefficient.
% DTM 29/1/2016

Fs = 44100;         % Sample Rate for audio output
T = 1/Fs;           % Sample period
N = 44100;          % Number of samples in output
f0 = 441;           % Fundamental Frequency of string in Hz.

L = floor(0.5*Fs/f0);   % String Length in samples - floor rounds each value to nearest int
x = 0.5;                % Pluck Position as a proportion of string length
pickup = floor(L/2);    % Output Position

r = -1;                 % Bridge Reflection Coefficient

% Right-going delay line, defined by L 
right = zeros(1,L);
% Left-going delay line, defined by L
left = zeros(1,L);

% Define initial string shape from pluck position x
pluck = x*(L-1); % Find point on string corresponding to pluck position
x = [ ([0:floor(pluck)]/pluck),(L-1-[(floor(pluck)+1):(L-1)])/(L-1-pluck)];

% Initial displacement for each delay line is equivalent to plucked string
% excitation shape divided by 2.
left(1:L) = x(1:L)/2;
right(1:L) = x(1:L)/2;

% Initialize output
out = zeros(1,N);

% Initialize variables for display
pkval = max(abs(x));
string_pos = 1:L;

% Main digital waveguide loop
for n = 1:N
    
    % Resamples the output at a rate of Fs/200 to plot the mass 
    % displacement in real-time.
    if mod(n,200)==0 
        % Plot left and right-moving waves, and their sum
        plot(string_pos, left, string_pos, right, string_pos, left+right);
        % Make sure the axis scaling stays the same for each plot
        axis([1 L -pkval pkval]);
        pause(0.02);
    end
    
    
  % Shift left-going wave one step left; append dummy value for now
  left = [left(2:L),0];
  % At the 'nut' (left-hand end), assume perfect reflection (* -1).
  % New right-going value is negative of new value at nut of left-going
  nut = -left(1);
  
  % Add some slight attenuation at the nut
  % nut = nut * 0.985;
  
  % Add reflection from nut into first element of right-going delay line;
  % Shift right-going wave one step
  right = [nut, right(1:L-1)];
  % At the 'bridge' (right-hand end), assume perfect reflection (* -1).
  % New left-going value is negative of new value at bridge of right-going
  % bridge = r*right(L);
  % bridge = r*0.5*(right(L) + right(L-1));
  
  % for an N point moving average filter:
  % bridge = r*(1/N)*(right(L) + right(L-1) + ... + right(L-(N-1));
   bridge = r*0.33*(right(L) + right(L-1) + right(L-2)); % 3 point moving average filter
  
  
  % Add a low pass filter to attenuate values at bridge
  % g = 0.5*right(L) + 0.5*right(L-1);
  % bridge = bridge * g;
  
  % Add new bridge value to end of left-going delay line, replacing dummy
  % value from above:
  left(L) = bridge;
  

  % Output is sum of left and right going delay lines at pickup point.
  % Calculate output:
  out(n) = left(pickup) + right(pickup);
 
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
axis([0 16000 -40 0]);
xlabel('Frequency (Hz)');
ylabel('Magnitude Response (dB)');

% Plot spectrogram of output
figure(4);
spectrogram(out, 'yaxis', 1024, 256, 1024);

soundsc(out, Fs);
