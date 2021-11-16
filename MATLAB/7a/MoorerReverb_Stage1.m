% The Moorer Reverb Design starting with the Feedback Comb Filter,
% developing into a Schroeder Reverb script to finally add
% Early Reflections and a low-pass filter in the feedback 
% path of the comb filters.
%
% Developed from examples in: 
% "Hack Audio: An Introduction to Computer Programming and Digital Signal
% Processing in MATLAB" by Eric Tarr.
%
% DTM, 6/11/2018
% 

clear;
clc;

% Load in an example test file
[in,Fs] = audioread('drums.wav');

% Zero pad the file so the output will include the reverb tail
in = [in;zeros(Fs*3,1)];

% Define impulse excitation for obtaining impulse response
IR_test = [zeros(Fs*3,1)];
IR_test(1) = 1;

% Initialize Main Output Signal
N_out = length(in);
out = zeros(N_out,1);

% Initialize Impulse Response Output Signal
N_IR = length(IR_test);
IR_out = zeros(N_IR,1);

% Set Maximum delay time for the unit reverberators of 70 ms
maxDelay = ceil(.07*Fs);

% Initialize all buffers - one buffer per unit reverberator
buffer1 = zeros(maxDelay,1); 

% Initialise the Early Reflection Unit Tapped Delay Line

% Delay (ms) and Gain Parameters
% Comb Filters
d1 = floor(.02*Fs); g1 = 0.75;

% Allpass filters

% Variables used as delay for a simple LPF in each Comb Filter function
fbLPF1 = 0;

% Impulse Response
for n = 1:N_IR
    
    % Early Reflection Tapped Delay Line
    
    % Parallel FBCFs
    [IR_out(n,1),buffer1,fbLPF1] = FeedbackComb(IR_test(n,1),buffer1,n,d1,g1,fbLPF1,false);
    
    % Two Series All-pass Filters
      
end
% Plot Impulse Response
figure();
plot(IR_out);

% Approximate Reverb Time:
% Only becomes useful with later examples.
IR_out = IR_out/max(abs(IR_out)); % Normalize to Unity Gain (0 dB)
Ts = 1/Fs;
t = [0:N_IR-1]*Ts;
figure();
plot(t,20*log10(abs(IR_out))); 
line([0 4],[-60 -60],'Color','red','LineStyle','--');
axis([0 4 -80 0]);

% Reverberate
for n = 1:N_out
    
    % Early Reflections Tapped Delay Line
    
    % Four Parallel FBCFs
    [out(n,1),buffer1,fbLPF1] = FeedbackComb(in(n,1),buffer1,n,d1,g1,fbLPF1,false);  
    
    % Two Series All-pass Filters
    
end
sound(out,Fs);



