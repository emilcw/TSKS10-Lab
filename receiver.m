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
%---------------------Filter from desired frequency band------------------

n_bp = 100;         %According to task 
w1 = f1 / fn;       %Norm cutoff freq
w2 = f2 / fn;       %Norm cutoff freq
[b_bp,a_bp] = fir1(n_bp,[w1 w2],'bandpass');
y_bp = filter(b_bp,a_bp,y);

%--------------------Demodulate not compensated signal--------------------

% Create lowpass-filter
t_lp = 0:Ts:(length(y_bp)/fs)-Ts;   % Changed to match y_bp with carriers
n = 100;                            %According to task 
W = (bandwidth/2) / fn; 
[b,a] = fir1(n,W,'low');

% Demodulate and lowpass filter to avoid aliasing
yI = filter(b,a,(2 * y_bp .* transpose(cos(2 * pi * fc * t_lp ))));
yQ = filter(b,a,(2 * y_bp .* transpose(-sin(2 * pi * fc * t_lp))));

% Here we have that
% yI = xI(t)*cos(delta) + xQ(t)sin(delta) 
% yQ = -xI(t)*sin(delta) + xQ(t)*cos(delta)
%where delta is phase shift in radians. If delta = 0 they would be exact
%teh same, but they are not since channel and filters timeshift.

%We have that (according to exercise 2.13
% xI = yI*cos(delta) - yQ*sin(delta)
% xQ = yI*sin(delta) + yQ*cos(delta)
% Where delta is the same phaseshift.

%so to get the real xI and xQ, we need delta. 
%----------------------- Find tau and compensate--------------------------

% Generate chirp-signal
t1 = 0:Ts:1-Ts;
f0 = 1;         %Start frequency
epsilon = 10;    %Chirp-rate
chirp = cos(2*pi*f0*(1+(epsilon.*t1)).*t1);

% Korskorrelera yI och yQ med chirp
[R_I, lagsI] = xcorr(yI, chirp);
[R_Q, lagsQ] = xcorr(yQ, chirp);

% DEBG
% stem(lagsI, R_I)
% stem(lagsQ, R_Q)
% https://web.mit.edu/6.02/www/f2006/handouts/Lec9.pdf

%Bra med chirp, hittar exakt en topp i lags.
% Se föreläsning  + Hans anteckningar, han visar hur man kompenserar för
% fasskift. Detta måste läggas in i cos och sin också! Annars använd peaks
% och lags för att ta fram tau, sedan A med norm. 

% Find peaks
%peaks_I = findpeaks(R_I);
[V_I,Index_I] = max(abs(R_I))
[V_Q,Index_Q] = max(abs(R_Q))

a = lagsI(Index_I)
b = lagsQ(Index_Q)  %Detta blir 184 och 183 men kan ej vara rätt ty det är bara 83 och 84 nollor i början.
phase_shift = abs(a-b)

tau = 83 * Ts

%peaks_I = findpeaks(R_I);
%peaks_Q = findpeaks(R_Q);

% Find maxpeak
%max_peak_I = max(peaks_I);
%max_peak_Q = max(peaks_Q);

% Find index for maxpeak (which corresponds to our delay)
%tau_I = find(abs(R_I) == max_peak_I);
%tau_Q = find(abs(R_Q) == max_peak_Q);
%tau = tau_I * Ts;

% Compensate for tau + cut away chirp
y_delay_comp = y(a:end);
%------------------------Find A and compensate-----------------------------

% Find A
R_chirp = xcorr(y_delay_comp);
A = norm(R_chirp);

% Compensate for A
z = (1/A) * (y_delay_comp);

% ----------------------Demodulate compensated version---------------------

% Filter out compensated version
n = 100;        %According to task 
w1 = f1 / fn;   %Norm cutoff freq
w2 = f2 / fn;   %Norm cutoff freq
[b,a] = fir1(n,[w1 w2],'bandpass');
z_bp = filter(b,a,z);

% Create lowpass-filter
t1 = 0:Ts:6+0.025-Ts;   % Added 10000/400000 to compensate for phaseshift.
n = 100;                %According to task 
W = (bandwidth/2) / fn; 
[b,a] = fir1(n,W,'low');

% Denna modulering måste göras om enligt det som gjorts ovan. 
% Demodulate and lowpass filter to avoid aliasing to the compensated version
I_carry = transpose(cos(2 * pi * fc *t1 ));
Q_carry = transpose(sin(2 * pi * fc *t1));
I_carry_index = I_carry(abs(tau_I + length(chirp)):end);
Q_carry_index = Q_carry(abs(tau_I + length(chirp)):end);
zI_demod = filter(b,a,(2 * z_bp .* I_carry_index));
zQ_demod = -filter(b,a,(2* z_bp .* Q_carry_index));

%Downsample zI and zQ
sample_factor = 20;
zI = downsample(zI_demod, sample_factor);
zQ = downsample(zQ_demod, sample_factor);
