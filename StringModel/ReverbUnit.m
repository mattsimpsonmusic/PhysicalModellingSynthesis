function [output] = ReverbUnit(input, room, mix)
% The Moorer Reverb Design starting with the Feedback Comb Filter,
% developing into a Schroeder Reverb script to finally add
% Early Reflections and a low-pass filter in the feedback 
% path of the comb filters.
%
% Developed from examples in: 
% "Hack Audio: An Introduction to Computer Programming and Digital Signal
% Processing in MATLAB" by Eric Tarr.
%
% Last edited on 31/12/2021
% Exam no: Y3858230


    % Function Arguments:
        % input                 input signal
        % room                  reverb enviroment preset
        % mix                   dry/wet mix

    % Define sample rate
    Fs = 44100;

    % Swap dimensions of input matrix to allow for effects processing
    in = input.';
    
    in = in / max(abs(in));   % Normalise to avoid clipping

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
    c1buffer = zeros(maxDelay,1);
    c2buffer = zeros(maxDelay,1);
    c3buffer = zeros(maxDelay,1);
    c4buffer = zeros(maxDelay,1);
    a1buffer = zeros(maxDelay,1);
    a2buffer = zeros(maxDelay,1);

    % Initialise the Early Reflection Unit Tapped Delay Line
    IR_comb_input = [zeros(Fs*3,1)];
    IR_comb_input(1) = 1;
    ERbuffer = zeros(maxDelay,1);
    comb_input = [in;zeros(Fs*3,1)];
    
    % Select delay and gain parameters based on preset %

    % Medium/large hall preset
    if room == 2
        % Delay (ms) and Gain Parameters
        % Comb Filters
        d1 = floor(.0297*Fs); g1 = 0.9;
        d2 = floor(.0371*Fs); g2 = 0.8;
        d3 = floor(.0411*Fs); g3 = 0.8;
        d4 = floor(.0437*Fs); g4 = 0.9;

        % Allpass filters
        a1 = floor(.005*Fs);  g5 = 0.9;
        a2 = floor(.0017*Fs); g6 = 0.8;
        
    % Very large chamber preset
    elseif room == 3
        % Delay (ms) and Gain Parameters
        % Comb Filters
        d1 = floor(9.79*Fs); g1 = 0.8;
        d2 = floor(6.71*Fs); g2 = 0.85;
        d3 = floor(5.11*Fs); g3 = 0.85;
        d4 = floor(4.37*Fs); g4 = 0.79;

        % Allpass filters
        a1 = floor(9*Fs);  g5 = 0.85;
        a2 = floor(4.7*Fs); g6 = 0.75;
        
    % Default preset (small room)
    else
        % Delay (ms) and Gain Parameters
        % Comb Filters
        d1 = floor(.0297*Fs); g1 = 0.50;
        d2 = floor(.0371*Fs); g2 = 0.65;
        d3 = floor(.0411*Fs); g3 = 0.45;
        d4 = floor(.0437*Fs); g4 = 0.50;

        % Allpass filters
        a1 = floor(.005*Fs);  g5 = 0.75;
        a2 = floor(.0017*Fs); g6 = 0.75;
    end

    % Variables used as delay for a simple LPF in each Comb Filter function
    fbLPF1 = 0.5;
    fbLPF2 = 0.5;
    fbLPF3 = 0.5;
    fbLPF4 = 0.5;

    % Impulse Response
    for n = 1:N_IR

        % Early Reflection Tapped Delay Line
        [IR_comb_input(n,1), ERbuffer] = EarlyReflections(IR_test(n,1), ERbuffer, Fs, n);

        % Parallel FBCFs
        [comb1,c1buffer,fbLPF1] = FeedbackComb(IR_comb_input(n,1),c1buffer,n,d1,g1,fbLPF1,true);
        [comb2,c2buffer,fbLPF2] = FeedbackComb(IR_comb_input(n,1),c2buffer,n,d2,g2,fbLPF2,true);
        [comb3,c3buffer,fbLPF3] = FeedbackComb(IR_comb_input(n,1),c3buffer,n,d3,g3,fbLPF3,true);
        [comb4,c4buffer,fbLPF4] = FeedbackComb(IR_comb_input(n,1),c4buffer,n,d4,g4,fbLPF4,true);

        % Sum parallel comb filters
        comb_IR_sum = 0.25*(comb1 + comb2 + comb3 + comb4);

        % Two Series All-pass Filters
        [a1_IR, a1buffer] = AllPass(comb_IR_sum, a1buffer, n, a1, g5);
        [a2_IR, a2buffer] = AllPass(a1_IR, a2buffer, n, a2, g6);

        % Assign output
        IR_out(n) = a2_IR;


    end

    % Plots %
    
    % Plot Impulse Response
    %figure();
    %plot(IR_out);

    % Approximate Reverb Time:
    %IR_out = IR_out/max(abs(IR_out)); % Normalize to Unity Gain (0 dB)
    %Ts = 1/Fs;
    %t = [0:N_IR-1]*Ts;
    %figure();
    %plot(t,20*log10(abs(IR_out))); 
    %line([0 4],[-60 -60],'Color','red','LineStyle','--');
    %axis([0 4 -80 0]);

    % Reverberate
    for n = 1:N_out

        % Early Reflections Tapped Delay Line
        [comb_input(n,1), ERbuffer] = EarlyReflections(in(n,1), ERbuffer, Fs, n);

        % Four Parallel FBCFs
        [comb1_out,c1buffer,fbLPF1] = FeedbackComb(comb_input(n,1),c1buffer,n,d1,g1,fbLPF1,true);
        [comb2_out,c2buffer,fbLPF2] = FeedbackComb(comb_input(n,1),c2buffer,n,d2,g2,fbLPF2,true);
        [comb3_out,c3buffer,fbLPF3] = FeedbackComb(comb_input(n,1),c3buffer,n,d3,g3,fbLPF3,true);
        [comb4_out,c4buffer,fbLPF4] = FeedbackComb(comb_input(n,1),c4buffer,n,d4,g4,fbLPF4,true);

        % Sum parallel comb filters
        comb_out_sum = 0.25*(comb1_out + comb2_out + comb3_out + comb4_out);

        % Two Series All-pass Filters
        [a1_out, a1buffer] = AllPass(comb_out_sum, a1buffer, n, a1, g5);
        [a2_out, a2buffer] = AllPass(a1_out, a2buffer, n, a2, g6);

        % Assign output
        out(n) = a2_out;

    end
    % Calculate output using dry/wet mix
    output = mix*out+(1-mix)*in;
end


