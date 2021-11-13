% Digital Waveguide Flute Example based on Perry Cook's Slide Flute
% Refed to Lab Scrip for the block diagram design
% Uses three-point moving average low pass reflection/loss filter 
% DTM 20/10/2020

Fs = 44100;                 % Sample Rate for audio output
N = 88200;                  % Number of samples in output
pressure = 0.5;             % Maximum Pressure input value
breath = 0.0085;            % Breath value
VibLevel = 0.05;            % Vibrato Amount
f0 = 440;                   % Fundamental Frequency in Hz.

attack = 0.025;             % attack time
sustain = 0.2;              % sustain start time 
release = 0.001;            % release time (Csound = 0.02)

as = round(N*attack);       % attack end sample
ss = round(N*sustain);      % sustain start = decay end sample
ds = ss-as;                 % decay time in samples
rs = round(N*(1-release));  % release start sample

Feedback1 = 0.7;            % Bore Delay to Embouchure Delay 
Feedback2 = 0.8;            % Bore Delay Reflection coefficient

boreDelay = 0.25*floor(Fs/f0);
JetDelay = floor(boreDelay/2);

% Flow Setup: ADSR envelope.
% attack
flow(1:as) = linspace(0,1.1*pressure,as); % was 1.1
% decay
flow(as+1:as+ds) = linspace(1.1*pressure,pressure,ds); %was 1.1
% sustain
flow(as+ds+1:rs) = pressure;
% release
flow(rs+1:N) = linspace(pressure,0,N-rs);

% Output Envelope
% attack
Env(1:1000) = linspace(0,1.0,1000);
% sustain
Env(1000+1:N-1000) = 1.0;
% release
Env(N-1000+1:N) = linspace(1.0,0,1000);

% Vibrato Envelope
VibEnv(1:N) = 1.0;

% Vibrato
vibrato = 0;

% Noise Excitation
noise = 2*rand(1,N)-1;

% Total Input
Excitation = (breath*(noise.*flow))+flow+vibrato*VibLevel.*VibEnv;

% Right-going delay line, defined by boreDelay 
right_bore = zeros(1,boreDelay);
% Left-going delay line, defined by boreDelay
left_bore = zeros(1,boreDelay);
%JetDelay
Jet = zeros(1,JetDelay);
% Initialize output
out = zeros(1,N);

% Main digital waveguide loop
for n = 1:N
      
    Jet(1) = Excitation(n) + Feedback1*left_bore(1);
    
    % Nonlinear interaction between excitation and flute bore
    r = Jet(JetDelay) - (Jet(JetDelay)*Jet(JetDelay)*Jet(JetDelay)); 
    if (r<-1) r=-1; end
    if (r>1) r=1; end

    right_bore(1) = r + Feedback2*left_bore(1);

    % Implement a three point moving average filter
    LPFilter = 0.33*(right_bore(boreDelay) + right_bore(boreDelay-1) + right_bore(boreDelay-2));
    left_bore(boreDelay) = -0.999*LPFilter;

    out(n) = right_bore(boreDelay);

    %Update Delay Lines
    Jet = [0, Jet(1:JetDelay-1)];
    right_bore = [0, right_bore(1:boreDelay-1)];
    left_bore = [left_bore(2:boreDelay),0];
    
end

%Apply Output Envelope to smooth input and output
out = out.*Env;

% Plot time domain response of output
figure(1);
plot(out);
xlabel('Time (samples)');
ylabel('Amplitude');

fftSize = 8192;
f = (0:fftSize-1)*(Fs/fftSize);

% Plot frequency domain response of output
figure(2);
semilogx(f,20*log10(abs(fft(out,fftSize))/max(abs(fft(out,fftSize)))));
axis([0 20000 -40 0]);
xlabel('Frequency (Hz)');
ylabel('Magnitude Response (dB)');

% Spectrogram to inspect effects of damping.
figure(3);
spectrogram(out,hann(1024),256,1024,Fs,'yaxis');

soundsc(out, Fs);
