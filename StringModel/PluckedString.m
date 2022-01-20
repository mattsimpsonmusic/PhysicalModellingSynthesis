function [out] = PluckedString(f0, x, pickup, r)
% Plucked String Digital Waveguide Sound Synthesis
% Based on travelling wave, two delay-line model
% Variable plucking position and output possible
% Simple bridge reflection coefficient.
% Last edited on 30/12/2021
% Exam no: Y3858230

    % Function Arguments:
    %   f0            fundamental frequency of the string in Hz.
    %   x             pluck position provide as a proportion of string length (useful range is 0<x<=1).
    %   pickup        output position (range is 0<pickup<=10) - recommended to use integers only.
    %   r             attenuation at the bridge termination of the string (range is 0<r<=1).

    % Initialise variables
    Fs = 44100;             % Sample Rate for audio output
    %T = 1/Fs;              % Sample period (not used)
    N = 44100;              % Number of samples in output
    L = floor(0.5*Fs/f0);   % String Length in samples - floor rounds each value to nearest int

    % Normalise pickup from input given (range of ints from 1-10)
    pickup = floor((pickup*L)/10);

    % Change sign of bridge attenuation coefficient for reflection
    r = -r;

    % Right-going delay line, defined by L 
    right = zeros(1,L);
    % Left-going delay line, defined by L
    left = zeros(1,L);

    % Define initial string shape from pluck position x
    pluck = x*(L-1); % Find point on string corresponding to pluck position
    x = [([0:floor(pluck)]/pluck),(L-1-[(floor(pluck)+1):(L-1)])/(L-1-pluck)];

    % Initial displacement for each delay line is equivalent to plucked string
    % excitation shape divided by 2.
    left(1:L) = x(1:L)/2;
    right(1:L) = x(1:L)/2;

    % Initialize output
    out = zeros(1,N);

    % Initialize variables for display (used for real time plotting)
    %pkval = max(abs(x));
    %string_pos = 1:L;

    % Main digital waveguide loop
    for n = 1:N
        
      % Shift left-going wave one step left; append dummy value for now
      left = [left(2:L),0];
      % At the 'nut' (left-hand end), assume perfect reflection (* -1).
      % New right-going value is negative of new value at nut of left-going
      nut = -left(1);

      % Add reflection from nut into first element of right-going delay line;
      % Shift right-going wave one step
      right = [nut, right(1:L-1)];
      % At the 'bridge' (right-hand end), assume perfect reflection (* -1).
      % New left-going value is negative of new value at bridge of right-going

      % Add a low pass filter to attenuate values at the bridge
      bridge = r*0.33*(right(L) + right(L-1) + right(L-2)); % 3-point moving average filter
      % bridge = r*right(L);                                % Simple reflection
      % bridge = r*0.5*(right(L) + right(L-1));             % 2-point moving average filter
      % for an N point moving average filter:
      % bridge = r*(1/N)*(right(L) + right(L-1) + ... + right(L-(N-1));

      % Add new bridge value to end of left-going delay line, replacing dummy
      % value from above:
      left(L) = bridge;

      % Output is sum of left and right going delay lines at pickup point:
      out(n) = left(pickup) + right(pickup);

    end
end

