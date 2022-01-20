function [output] = AmpUnit(input, distortion_gain, drive)
% A Simple Guitar Amp Model
% This takes the distortion and overdrive effects from the DAFx book
% and combines them with the impulse response of a miked up guitar amp
% 4x12 speaker cabinet to give a more than passable guitar amp model.
% Last edited on 30/12/2021
% Exam no: Y3858230

    % Function Arguments:
        % input                 input signal
        % distortion_gain       gain for distortion unit
        % drive                 overdrive

    % Swap dimensions of input matrix to allow for effects processing
    audio_input = input.'; 

    % Distortion Block
    out_dist = distortion(audio_input,distortion_gain,1);

    % Overdrive Block
    out_drive = overdrive(drive*out_dist);

    % Cabinet Section
    [CabIR, ~] = audioread('Cab412_IR.wav'); % Import Cabinet simulation IR
    out_cab = conv(out_drive, CabIR);

    output = out_cab;
end

