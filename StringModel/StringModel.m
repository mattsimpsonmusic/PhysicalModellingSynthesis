function [] = StringModel(f0, x, pickup, r, type, reverse, amp, gain, drive, reverb, room, mix, chorus, wobble, feedback, membrane)
% A Functional String Model Designed for Sound Synthesis and Sound Design Purposes
% This is the top level function to run the physical model of the
% synthesised string.
%
% Last edited on 09/01/2022
% Exam no: Y3858230

    % Function Arguments %
        % f0            (fundamental frequency)
        % x             (pluck position - plucked string only)
        % pickup        (pickup position - plucked string only)
        % r             (damping - plucked string only)
        % type          (string type & chord arrangement)
        % reverse       (reverse output sound)
        % amp           (amp unit on/off)
        % gain          (amp distortion gain)
        % drive         (amp overdrive)
        % reverb        (reverb unit on/off)
        % room          (reverb environment)
        % mix           (reverb dry/wet mix)
        % chorus        (chorus unit on/off)
        % wobble        (chorus modulation rate)
        % feedback      (chorus feedback percentage)
        % membrane      (membrane excitation on/off)
        
        
    % Function Argument Descriptions %
    
    % String Parameters %
    
        % 'f0' is the fundamental frequency of the string in Hz. 

        % 'x' is the pluck position provided as a proportion of string length 
        % (useful range is 0<x<=1).

        % 'pickup' is the output position of the plucked string model.
        % (range is 0<pickup<=10) - recommended to use integers only.

        % 'r' is the attenuation coefficient at the bridge termination of the string 
        % (range is 0<r<=1).

        % 'type' selects between the different string types (odd numbers select
        % plucked string model and even select a struck string piano model),
        % with different argument numbers selecting a different type of sound
        % (options are a single note, major chord, minor chord, and diminished chord).
    
    % Effects units %
    
        % 'reverse' provides the option to reverse the sound generated from
        % the string model. A value of '1' will reverse the string sound.
        
        % 'amp' controls the on/off switch for the amplification unit -
        % selecting a '1' will pass the string through the amplification
        % unit.
        
        % 'gain' provides the amount of distortion gain inside the
        % amplification unit. Values above '0' will work, but a suitable
        % range for most cases would be between '1' and '15'.
        
        % 'drive' is the amount of overdrive present inside the
        % amplification unit. Values above '0' will work, but a suitable
        % range for most cases would be between '1' and '15'.
        
        % 'reverb' controls the on/off switch for the reverberation unit -
        % selecting a '1' will pass the string through the reverberation
        % unit.
        
        % 'room' is a parameter which provides the option of three
        % different pre-calculated simulated enviromnents for the user -
        % '1' will choose a small room, '2' a medium/large hall, and '3'
        % for a very large chamber/cathedral enviroment.
        
        % 'mix' changes the dry/wet mix of the reverberation unit using
        % values between 0 and 1, where 0 is the original signal (dry) and
        % 1 is the fully reverberated signal (wet).
        
        % 'membrane' controls whether or not the string is loaded and
        % convolved with a 2D square membrane. A value of '1' will load the
        % string onto the membrane.
        
        % 'chorus' determines whether or not choruse is applied to the
        % string sound - a value of '1' will apply chorus to the string.
        
        % 'wobble' sets the modulation frequency of the chorus unit in Hz.
        % This will take any value above '0', but using larger numbers will
        % produce some interesting transofrmations of the output sound.
        
        % 'feedback' controls the percentage of feedback from the chorus
        % unit in the signal. A value of '0.0' will have no feedback and a
        % value of '1.0' will set the maximum feedback value.
        
        % 'membrane' will toggle on and off the application of using the
        % string as an excitation to a pre-defined 2D membrane, and then
        % convolving the output of that membrane with the original string
        % signal. A value '1' will set the unit on.


    % Default Function Settings % 
    
        % The command outlined below will play a plucked string at A4 with no
        % additional effects, heard half way along the string:

        % StringModel(440, 0.5, 5, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        

    %%%% Main Function %%%%

    Fs = 44100; % Sample rate for audio playback
   

    % Selection process to choose between string type
    % Single plucked string
    if type == 1
        string = PluckedString(f0,x,pickup,r);
    % Single struck string
    elseif type == 2
        string = Piano(f0);
    % Plucked string major chord
    elseif type == 3
        string1 = PluckedString(f0,x,pickup,r);         % Tonic
        string2 = PluckedString(f0*1.256,x,pickup,r);   % Maj3
        string3 = PluckedString(f0*1.498,x,pickup,r);   % Fifth
        string = (string1 + string2 + string3)/3;
    % Struck string major chord
    elseif type == 4
        string1 = Piano(f0);                            % Tonic
        string2 = Piano(f0*1.256);                      % Maj3
        string3 = Piano(f0*1.498);                      % Fifth
        string = (string1 + string2 + string3)/3;
    % Plucked string minor chord
    elseif type == 5
        string1 = PluckedString(f0,x,pickup,r);         % Tonic
        string2 = PluckedString(f0*1.189,x,pickup,r);   % Min3
        string3 = PluckedString(f0*1.498,x,pickup,r);   % Fifth
        string = (string1 + string2 + string3)/3;
    % Struck string minor chord
    elseif type == 6
        string1 = Piano(f0);                            % Tonic
        string2 = Piano(f0*1.189);                      % Min3
        string3 = Piano(f0*1.498);                      % Fifth
        string = (string1 + string2 + string3)/3;
    % Plucked string diminished chord
    elseif type == 7
        string1 = PluckedString(f0,x,pickup,r);         % Tonic
        string2 = PluckedString(f0*1.189,x,pickup,r);   % Min3
        string3 = PluckedString(f0*1.414,x,pickup,r);   % Sharp4
        string4 = PluckedString(f0*1.682,x,pickup,r);   % Sixth
        string = (string1 + string2 + string3 + string4)/4;
    % Struck string diminished chord
    elseif type == 8
        string1 = Piano(f0);                            % Tonic
        string2 = Piano(f0*1.189);                      % Min3
        string3 = Piano(f0*1.414);                      % Sharp4
        string4 = Piano(f0*1.682);                      % Sixth
        string = (string1 + string2 + string3 + string4)/4;
    else 
        error("Invalid type entered - please select a number between 1-8 to set string type");
    end


    % Reverse playback section
    if reverse == 1
        string = fliplr(string); % Flips every element in the array back to front
    end

    % Amplification section
    if amp == 1
        string = AmpUnit(string, gain, drive);
    end

    % Membrane Section
    if membrane == 1
        membraneout = Membrane(string);
        string = conv(string, membraneout);
    end
    
    % Reverb section
    if reverb == 1
        string = ReverbUnit(string, room, mix);
    end

    % Chorus Section
    if chorus == 1
        string = ChorusUnit(string, wobble, feedback);
    end

    % Plot String Output %

    % Plot time domain response of output
    figure(1);
    plot(string);
    xlabel('Time (samples)');
    ylabel('Amplitude');

    % Plot frequency domain response of output
    fftSize = 8192;
    f = (0:fftSize-1)*(Fs/fftSize);
    figure(2);
    semilogx(f,20*log10(abs(fft(string,fftSize))/max(abs(fft(string,fftSize)))));
    axis([0 16000 -40 0]);
    % amp axis - axis([0 20000 -100 0]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude Response (dB)');
    

    % Plot spectrogram of output
    figure(3);
    spectrogram(string, 'yaxis', 1024, 256, 1024);

    
    % Play String Output % 
    soundsc(string, Fs);
    
    % Write Output to Audio File %
    stringoutput = string / max(abs(string));   % Normalise to avoid clipping
    audiowrite('StringModel.wav', stringoutput, Fs); % Write to WAVE file
end

