function [output] = Membrane(input)
% 2D membrane Absorbing
% Note Courant Lambda value = sqrt(0.5)
% DTM (10/2020)
%
% Last edited on 01/01/2022
% Exam no: Y3858230

    % Function Arguments:
        % input                 input signal

    % MESH DEFINITION %

    sr = 48000;     % sample rate of system/mesh

    Nx = 100;       % number of grid points in x direction
    Ny = 100;       % number of grid points in y direction
    Ns = 4000;      % number of samples

    Inx = 50;       % Define an excitation point
    Iny = 50;       

    Outx = 10;      % Define an output point
    Outy = 10;

    % CALCULATE SIMPLIFIED CONSTANTS A,B & C

    lambda = sqrt(0.5);

    A = 2 - 4*(lambda^2);
    B = lambda^2;
    C = -1;

    % INITIALISE VARIABLES

    % Initialise matrices used to store values at t+1, t, and t-1
    p = zeros(Nx,Ny);       % p(t+1) at new time instant
    p1 = zeros(Nx,Ny);      % p(t) at current time instant
    p2 = zeros(Nx,Ny);      % p(t-1) at previous time instant

    % Initialise an array to capture output at a single point
    out = zeros(1, Ns);

    % SET EXCITATION

    % First define a 9x9 point region in the centre of the mesh
    Ne = 4;
    xe = (Inx-Ne):(Inx+Ne);
    ye = (Iny-Ne):(Iny+Ne);
    % Transform input into a 2D function/shape
    shape = input*input';
    % Load the shape into the centre of the mesh at current timestep (t) 
    p1(xe,ye) = shape;
    % Let the mesh at the previous timestep, t-1 have the same shape.
    p2 = p1;

    % MESH ANIMATION PLOT %

    % Draw the mesh at this starting point, t=0
    %figure(1);
    %clf;
    %surf(p1);
    %axis off;
    %axis equal;
    %shading interp;
    %colormap(copper);
    %V =  [1 100 1 100 -25 25];
    %axis(V);
    %pause; %You need to hit a key on the keyboard to start the animation.

    g = 0.01;    % This quantity that controls the boundary absorption 
    rt2 = sqrt(2.0);

    % Coefficients used in the absorbing boundary equations

    a0 = rt2/(rt2 + g);
    a1 = 1/(2.0 + (rt2*g));
    a2 = (g - rt2)/(g+rt2);

    % MAIN LOOP
    for n=1:Ns
        % finite difference equation

        for l = 1:Nx                              
            for m = 1:Ny   
                if (l == 1)&&(m == 1)          % Bottom Left Corner  
                    p(l,m) = 0;
                elseif (l == Nx)&&(m == 1)     % Bottom Right Corner
                    p(l,m) = 0;
                elseif (l == 1)&&(m == Ny)     % Bottom Right Corner
                    p(l,m) = 0;
                elseif (l == Nx)&&(m == Ny)     % Bottom Right Corner
                    p(l,m) = 0;
                elseif (l == 1)&&(m~=Ny)&&(m~=1)             % Left Boundary
                    p(l,m) = a0*p1(l+1,m) + a1*(p1(l,m+1) + p1(l,m-1))+ a2*p2(l,m);
                elseif (l == Nx)&&(m~=Ny)&&(m~=1)            % Right Boundary
                    p(l,m) = a0*p1(l-1,m) + a1*(p1(l,m+1) + p1(l,m-1))+ a2*p2(l,m);
                elseif (m == 1)&&(l~=Nx)&&(l~=1)             % Bottom Boundary
                    p(l,m) = a0*p1(l,m+1) + a1*(p1(l+1,m) + p1(l-1,m))+ a2*p2(l,m);
                elseif (m == Ny)&&(l~=Nx)&&(l~=1)            % Top Boundary
                    p(l,m) = a0*p1(l,m-1) + a1*(p1(l+1,m) + p1(l-1,m))+ a2*p2(l,m);
                elseif (m~=1)&&(m~=Ny)&&(l~=1)&&(l~=Nx) 
                     % main finite difference expression
                     p(l,m) = (A * p1(l,m)) + (B * (p1(l+1, m) + p1(l-1, m) + p1(l, m+1) + p1(l, m-1))) + (C * p2(l,m));
                end
            end
        end
        out(n) = p(Outx,Outy);
        % update mesh history
        p2 = p1;
        p1 = p;

        % MESH ANIMATION PLOT UPDATES %
        
        % update plot
        %if mod(n,50)==0 
        %    surf(p);
        %    axis off;
        %    axis equal;
        %    axis(V);
        %    shading interp;
        %    pause(0.001);
        %end
    end

    % Plot time domain response of output point
    figure(2);
    clf;
    plot(out);
    xlabel('time (samples)');
    ylabel('Amplitude');

    fftSize = 8192;
    f = (0:fftSize-1)*(sr/fftSize);

    % Plot frequency domain response of output point
    figure(3);
    clf;
    semilogx(f,20*log10(abs(fft(out,fftSize))/max(abs(fft(out,fftSize)))));
    axis([0 8000 -40 0]);
    xlabel('frequency (Hz)');
    ylabel('magnitude response (dB)');

    h = fir1(20, 0.25);  % 20th order low-pass filter with cutoff at 0.25* the
                         % sample rate normalised to half-H=Nyquist  
    % Create filtered output                     
    output = filter(h,1,out);  % filtered version of 'out' using 'h'

end
