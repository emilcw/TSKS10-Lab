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

%Compensate for filter delay
y_bp = y_bp((n_bp/2)+1:end);

%--------------------Demodulate not compensated signal--------------------

% Create lowpass-filter
t_lp = 0:Ts:(length(y_bp)/fs)-Ts;   % Changed to match y_bp with carriers
n_lp = 100;                            % According to task 
W = (bandwidth/2) / fn; 
[b_lp,a_lp] = fir1(n_lp, W,'low');

%Create Carriers
I_carry = transpose(cos(2 * pi * fc * t_lp ));
Q_carry = transpose(sin(2 * pi * fc * t_lp ));

% Create functions to be demodulated
yI_before = 2 * y_bp .* I_carry;
yQ_before = 2 * y_bp .* Q_carry;

% Demodulate and lowpass filter to avoid aliasing
yI = filter(b_lp,a_lp, yI_before);
yQ = -filter(b_lp,a_lp,yQ_before);

% Compensate for filter delay
yI = yI((n_lp/2):end);
yQ = yQ((n_lp/2):end);

%----------------------- Find tau and compensate--------------------------

% Generate chirp-signal
t1 = 0:Ts:1-Ts;
f0 = 1;         %Start frequency
epsilon = 10;    %Chirp-rate
chirp = cos(2*pi*f0*(1+(epsilon.*t1)).*t1);

% Korskorrelera yI och yQ med chirp
[R_I, lagsI] = xcorr(yI, chirp);
[R_Q, lagsQ] = xcorr(yQ, chirp);

% Find peaks
[V_I,Index_I] = max(abs(R_I));
[V_Q,Index_Q] = max(abs(R_Q));

% Find tau
tau_I = lagsI(Index_I)
tau_Q = lagsQ(Index_Q)
tau = tau_I * Ts;

figure(5)
plot(yI)
title('not comp')

% Compensate for tau
y_comp_tau = y_bp((tau_I+1):end); % y where tau is compensated

figure(6)
plot(y_comp_tau)
title('comp')

%------------------------Find A and compensate-----------------------------

% Find A
R_C = xcorr(chirp, y_comp_tau); 
cor_top = norm(R_C);            % top of correlation
chirp_energy = norm(chirp)^2;
A = cor_top / chirp_energy;

% Compensate for A
z = (1/A) * (y_comp_tau);

% Cut away chirp
%z = z(length(chirp)+1:end);

% ----------------------Demodulate compensated version---------------------

% Filter out compensated version
%n_bp = 100;        %According to task 
%w1 = f1 / fn;   %Norm cutoff freq
%w2 = f2 / fn;   %Norm cutoff freq
%[b,a] = fir1(n,[w1 w2],'bandpass');
%z_bp = filter(b,a,z);

% Compensate for filter delay
%z_bp = z_bp((n_bp/2)+1:end);

% Create lowpass-filter
t1 = 0:Ts:(length(z)/fs)-Ts;
n = 100;                %According to task 
W = (bandwidth/2) / fn; 
[b,a] = fir1(n,W,'low');

% Demodulate and lowpass filter to avoid aliasing to the compensated version

%Create carries
I_carry = transpose(cos(2 * pi * fc *t1 ));
Q_carry = transpose(sin(2 * pi * fc *t1));

% LP filter before downsample
zI_demod = filter(b,a,(2 * z .* I_carry));
zQ_demod = -filter(b,a,(2* z .* Q_carry));

%Compensate for filter delay
zI_demod = zI_demod((n_lp/2):end);
zQ_demod = zQ_demod((n_lp/2):end);

%Get the real values
%delta = tau;
%zI_demod = zI_demod*cos(delta) + zQ_demod*sin(delta);
%zQ_demod = zI_demod*sin(delta) + zQ_demod*cos(delta);

%Downsample zI and zQ
sample_factor = 20;
zI = downsample(zI_demod, sample_factor);
zQ = downsample(zQ_demod, sample_factor);

zI = zI(1:100000);
zQ = zQ(1:100000);
end
