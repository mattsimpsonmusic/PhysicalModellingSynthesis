function [output] = ChorusUnit(input, modulation_rate, feedback)
% A simple chorus effects unit, originally written by Dan Prince
% Orignal source:
% https://github.com/danpprince/matlab-chorus
% Last edited on 31/12/2021
% Exam no: Y3858230

    % Function Arguments:
        % input                 input signal
        % modulation_rate       frequency of modulation in Hz
        % feedback              percentage of feedback in output signal
        
    % Define sample rate
    Fs = 44100;

    % Initialise variables
    delay_length     = 0.013;   % delay of chorus effect (seconds)
    modulation_depth = 0.003;   % depth of signal modulation (seconds)
    low_shelf_freq   = 600;     % frequency of shelf (Hz)
    low_shelf_gain   = -7;      % shelf gain (dB)
    mix = 1;                    % dry/wet mix (set to full wet)

    % Calculate delay time and modulation depth
    delay_length_samples     = round(delay_length * Fs);
    modulation_depth_samples = round(modulation_depth * Fs);

    % Create delay buffer and output buffer
    modulated_output = zeros(length(input), 1);
    delay_buffer     = zeros(delay_length_samples + modulation_depth_samples, 1);

    % Argument for sin() modulation function. Converts the loop's control variable into 
    % the appropriate argument in radians to achieve the specified modulation rate
    modulation_argument = 2 * pi * modulation_rate / Fs;

    % Apply modulation
    for i = 1:(length(input))
        % Find index to read from for modulated output
        modulated_sample = modulation_depth_samples * sin(modulation_argument * i);
        modulated_sample = modulated_sample + delay_length_samples;

        % Get values to interpolate between
        interp_y1 = delay_buffer(floor(modulated_sample));
        interp_y2 = delay_buffer( ceil(modulated_sample));

        % Find difference between sample value and nearest integer
        query_sample = modulated_sample - floor(modulated_sample);

        % Interpolate to find the output value
        modulated_output(i) = interp_y1 + (interp_y2 - interp_y1) * (query_sample);

        % Save the input's current value in the buffer and advance to the next value
        new_sample = (input(i) + modulated_output(i) * feedback);
        delay_buffer = [ new_sample; delay_buffer(1 : length(delay_buffer)-1) ];
    end

    % Create low shelf filter
    % originally from http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
    w0     = 2 * pi * low_shelf_freq / Fs;
    S      = 0.5;
    A      = 10 ^ (low_shelf_gain / 40);
    alpha  = sin(w0) / 2 * sqrt( (A + 1/A) * (1/S - 1) + 2 );

    b0 =    A*( (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha );
    b1 =  2*A*( (A-1) - (A+1)*cos(w0)                   );
    b2 =    A*( (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha );
    a0 =        (A+1) + (A-1)*cos(w0) + 2*sqrt(A)*alpha;
    a1 =   -2*( (A-1) + (A+1)*cos(w0)                   );
    a2 =        (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha;


    % Apply low shelf EQ to the modulated signal
    modulated_output = filter([b0, b1, b2], [a0, a1, a2], modulated_output);

    % Add the dry and wet signals to get the final mixed version
    output = ((1 - mix) * input(:, 1) ) + (mix * modulated_output);

end

