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
buffer2 = zeros(maxDelay,1);
buffer3 = zeros(maxDelay,1);
buffer4 = zeros(maxDelay,1);
buffer5 = zeros(maxDelay,1);
buffer6 = zeros(maxDelay,1);

% Initialise the Early Reflection Unit Tapped Delay Line

% Delay (ms) and Gain Parameters
% Comb Filters
d1 = floor(.0297*Fs); g1 = 0.75;
d2 = floor(.0371*Fs); g2 = 0.75;
d3 = floor(.0411*Fs); g3 = 0.75;
d4 = floor(.0437*Fs); g4 = 0.75;

% Allpass filters
a1 = floor(.005*Fs);  g5 = 0.75;
a2 = floor(.0017*Fs); g6 = 0.75;

% Variables used as delay for a simple LPF in each Comb Filter function
fbLPF1 = 0;
fbLPF2 = 0;
fbLPF3 = 0;
fbLPF4 = 0;

% Impulse Response
for n = 1:N_IR
    
    % Early Reflection Tapped Delay Line
    
    % Parallel FBCFs
    [comb1,buffer1,fbLPF1] = FeedbackComb(IR_test(n,1),buffer1,n,d1,g1,fbLPF1,false);
    [comb2,buffer2,fbLPF2] = FeedbackComb(IR_test(n,1),buffer2,n,d2,g2,fbLPF2,false);
    [comb3,buffer3,fbLPF3] = FeedbackComb(IR_test(n,1),buffer3,n,d3,g3,fbLPF3,false);
    [comb4,buffer4,fbLPF4] = FeedbackComb(IR_test(n,1),buffer4,n,d4,g4,fbLPF4,false);
    
    % Sum parallel comb filters
    comb_IR_sum = 0.25*(comb1 + comb2 + comb3 + comb4);
      
    % Two Series All-pass Filters
    [a1_IR, buffer5] = AllPass(comb_IR_sum, buffer5, n, a1, g5);
    [a2_IR, buffer6] = AllPass(a1_IR, buffer6, n, a2, g6);
    
    % Assign output
    IR_out(n) = a2_IR;
    
      
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
    [comb1_out,buffer1,fbLPF1] = FeedbackComb(in(n,1),buffer1,n,d1,g1,fbLPF1,false);
    [comb2_out,buffer2,fbLPF2] = FeedbackComb(in(n,1),buffer2,n,d2,g2,fbLPF2,false);
    [comb3_out,buffer3,fbLPF3] = FeedbackComb(in(n,1),buffer3,n,d3,g3,fbLPF3,false);
    [comb4_out,buffer4,fbLPF4] = FeedbackComb(in(n,1),buffer4,n,d4,g4,fbLPF4,false);
    
    % Sum parallel comb filters
    comb_out_sum = 0.25*(comb1_out + comb2_out + comb3_out + comb4_out);
    
    % Two Series All-pass Filters
    [a1_out, buffer5] = AllPass(comb_out_sum, buffer5, n, a1, g5);
    [a2_out, buffer6] = AllPass(a1_out, buffer6, n, a2, g6);
    
    % Assign output
    out(n) = a2_out;
    
end

audiowrite('drums_7c.wav', out, Fs);
sound(out,Fs);




