% 2D membrane clamped
% Note Courant Lambda value = sqrt(0.5)
% DTM (10/2020)

% MESH DEFINITION

sr = 48000;     % sample rate of system/mesh

Nx = 100;       % number of grid points in x direction
Ny = 100;       % number of grid points in y direction

Ns = 4000;      % number of samples

Inx = 50;       % Define an excitation point
Iny = 50;       

Outx = 40;      % Define an output point
Outy = 40;

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
% Define a hanning window function, with an amplitude of 5.0 and a window
% size of 9 sample points - this is essenitally a smoothed impulse.
sh =  5*hanning(Ne*2+1);
% Transform this into a 2D function/shape
shape = sh*sh';
% Load the shape into the centre of the mesh at current timestep (t) 
p1(xe,ye) = shape;
% p1(xe, ye) = 343/2 * sqrt(1+1); % half cycle sine wave
% p1(xe+10, ye-5) = shape;  % add extra plucks
% p1(xe-20, ye-26) = shape;
% Let the mesh at the previous timestep, t-1 have the same shape.
p2 = p1;

% PLUCK %

% % First define a 9x9 point region in the centre of the mesh
% Ne = 4;
% xe = (Inx-Ne):(Inx+Ne);
% ye = (Iny-Ne):(Iny+Ne);
% % Define a hanning window function, with an amplitude of 5.0 and a window
% % size of 9 sample points - this is essenitally a smoothed impulse.
% sh =  5*hanning(Ne*2+1);
% % Transform this into a 2D function/shape
% shape = sh*sh';
% % Load the shape into the centre of the mesh at current timestep (t) 
% p1(xe,ye) = shape;
% % p1(xe, ye) = 343/2 * sqrt(1+1); % half cycle sine wave
% % p1(xe+10, ye-5) = shape;  % add extra plucks
% % p1(xe-20, ye-26) = shape;
% % Let the mesh at the previous timestep, t-1 have the same shape.
% p2 = p1;


% SINE WAVE %

% First define a region containing the entirity of thhe mesh
xe = 1:Nx;
ye = 1:Ny;
% Define a half cycle sine wave corresponding to the lowest mode of
% vibration
input_wave_x = 5*sin((pi*xe)/size(xe, 2));
input_wave_y = 5*sin((pi*ye)/size(ye, 2));
% Transform this into a 2D function/shape
shape = input_wave_x.*input_wave_y';
% Load the shape into the centre of the mesh at current timestep (t) 
p1(xe,ye) = shape;
% Let the mesh at the previous timestep, t-1 have the same shape.
p2 = p1;




% MESH ANIMATION

% Draw the mesh at this starting point, t=0
figure(1);
clf;
surf(p1);
axis off;
axis equal;
shading interp;
colormap(copper);
V =  [1 100 1 100 -25 25];
axis(V);
pause();  %You need to hit a key on the keyboard to start the animation.

% MAIN LOOP
for n=1:Ns
    % finite difference equation
    for l = 2:Nx-1                              
        for m = 2:Ny-1   
            p(l,m) = (A * p1(l,m)) + (B * (p1(l+1, m) + p1(l-1, m) + p1(l, m+1) + p1(l, m-1))) + (C * p2(l,m));
        end
    end
    
    % Get Output Value
    
    out(n) = p(Outx,Outy);
   
    % Update mesh history
    
    p2 = p1;
    p1 = p;
    
    % update plot
    if mod(n,200)==0 
        surf(p);
        axis off;
        axis equal;
        axis(V);
        shading interp;
        pause(0.01);
    end
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

% Play the output
h = fir1(20, 0.25);  % 20th order low-pass filter with cutoff at 0.25* the
                     % sample rate normalised to half-H=Nyquist                
lowpassout = filter(h,1,out);  % filtered version of 'out' using 'h'
sound(lowpassout, sr, 16);  % plays vector as a sound
