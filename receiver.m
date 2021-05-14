function [zI,zQ,A,tau] = receiver(y)
% RECEIVER        Given a signal y from a channel we demodulate it
%                 to find its initial message.
%
%   RECEIVER(y) Given a signal y(t) from a channel 
%               we filter it out from a desired frequency 
%               band, determines its amplitude and delay, 
%               compensates for it and demodulates 
%               its compensated version.

%------------------------Internal variables-------------------------------
f1 = 85e3;
f2 = 105e3;
fc = (f1 + f2) / 2;
bandwidth = 5e3;
fs = 400e3;
fn = fs / 2;
Ts = 1/fs;
signal_length = 100000;

%---------------------Filter from desired frequency band------------------

n_bp = 100;         %According to task 
w1 = f1 / fn;       %Norm cutoff freq
w2 = f2 / fn;       %Norm cutoff freq
[b_bp,a_bp] = fir1(n_bp,[w1 w2],'bandpass');
y_bp = filter(b_bp,a_bp,y);

% Compensate for filter delay
y_bp = y_bp((n_bp/2)+1:end);

%--------------------Demodulate not compensated signal--------------------

% Create Carriers
t1 = 0:Ts:(length(y_bp)*Ts)-Ts;  
I_carry = transpose(cos(2 * pi * fc * t1 ));
Q_carry = transpose(sin(2 * pi * fc * t1 ));

% Create functions to be demodulated
yI_demod = 2 * y_bp .* I_carry;
yQ_demod = 2 * y_bp .* Q_carry;

% Create lowpass-filter
n_low = 100; % According to task 
W = (2 * bandwidth)/ fn; 
[b_low,a_low] = fir1(n_low, W,'low');

% Demodulate and lowpass filter to avoid aliasing
yI = filter(b_low,a_low, yI_demod);
yQ = -filter(b_low,a_low,yQ_demod);

% Compensate for filter delay
yI = yI((n_low/2)+1:end);
yQ = yQ((n_low/2)+1:end);

%----------------------- Find tau and compensate--------------------------

% Generate chirp-signal
t2 = 0:Ts:1-Ts;
f0 = 1;          %Start frequency
epsilon = 10;    %Chirp-rate
chirp = cos(2*pi*f0*(1+(epsilon.*t2)).*t2);

% Add zeros to chirp to fix interference from channel
zeros_to_add = zeros(1, 40);
chirp_zeros = [chirp zeros_to_add];

% Cross-correlation with chirp to find tau
[R_I, lagsI] = xcorr(yI, chirp_zeros);
[R_Q, lagsQ] = xcorr(yQ, chirp_zeros);

% Find peakvalue and index
[V_I,Index_I] = max(abs(R_I));
[V_Q,Index_Q] = max(abs(R_Q));

% Determine the right tau, in other words the highest peak.
if V_I >= V_Q
    tau_index = lagsI(Index_I);
else
    tau_index = lagsQ(Index_Q);
end

% Convert tau_index to time in seconds.
tau = tau_index * Ts;

% Compensate for tau
z = y_bp(tau_index+1:end);

% ----------------Demodulate version where tau is compensated--------------

% Create lowpass-filter
t2 = 0:Ts:(length(z)*Ts)-Ts;
n_low2 = 100;                %According to task 
W = (2 * bandwidth) / fn; 
[b_low2,a_low2] = fir1(n_low2,W,'low');

% Create carries
I_carry = transpose(cos(2 * pi * fc * t2 ));
Q_carry = transpose(sin(2 * pi * fc * t2));

% Demodulate and...
zI_demod = 2 * z .* I_carry;
zQ_demod = 2 * z .* Q_carry;

% Lowpass-filter to avoid aliasing (and before downsampling)
zI_low = filter(b_low2,a_low2,zI_demod);
zQ_low = -filter(b_low2,a_low2,zQ_demod);

% Compensate for filter delay
zI_low = zI_low((n_low2/2)+1:end);
zQ_low = zQ_low((n_low2/2)+1:end);

%--------------Find A, compensate and cut away chirp-----------------------

% Find A
R_C = xcorr(chirp_zeros, zI_low); 
correlation_max = max(abs(R_C));           
chirp_energy = norm(chirp_zeros)^2;
A = correlation_max / chirp_energy;

% Cut off chirp
zI_no_chirp = zI_low(length(chirp_zeros)+1:end);
zQ_no_zeros = zQ_low(length(chirp_zeros)+1:end);

% Compensate for A
zI_no_chirp = (1/A) * (zI_no_chirp);
zQ_no_zeros = (1/A) * (zQ_no_zeros);

% Downsample zI and zQ (Lowpass filtering done above to interpolate)
sample_factor = 20;
zI_downsampled = downsample(zI_no_chirp, sample_factor);
zQ_downsampled = downsample(zQ_no_zeros, sample_factor);

% Compensate for signals beging mixed togheter. 
delta = max(angle(zQ_downsampled./zI_downsampled));
delta2 = mod(delta, 2*pi);
zI_fixed = zI_downsampled * cos(delta2) - zQ_downsampled * sin(delta2);
zQ_fixed = zQ_downsampled * sin(delta2) + zQ_downsampled * cos(delta2);
    
% Cut out the signal to its original length
zI = zI_fixed(1:signal_length);
zQ = zQ_fixed(1:signal_length);

end
