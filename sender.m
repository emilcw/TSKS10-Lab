function x = sender(xI, xQ)
% SENDER     Upsample xI, xQ, apply chirp and generate 
%            signal ready to be sent to channel.
%
%   SENDER(xI, xQ) returns x(t) upsampled with an applied 
%                  pulseform and the beginning
%                  so that it can be sent over a channel 
%                  via I/Q- modulation.
%    
%-------------------------Internal variables-------------------------------
f1 = 85e3;
f2 = 105e3;
bandwidth = 5e3;
fs = 400e3;
fn = fs / 2;
Ts = 1/fs;
fc = (f1 + f2) / 2;
sample_factor = 20;

%------------------------ Upsample to 400kHz-------------------------------

% Multiply by sample_factor to compensate for loss of signal energy
xI_upsamp = upsample(xI, sample_factor) * sample_factor;
xQ_upsamp = upsample(xQ, sample_factor) * sample_factor;

% Lowpass-filter to make it interpolate
n_low = 100;                %According to task 
W = (2 * bandwidth) / (fn);
[b_low,a_low] = fir1(n_low,W,'low');
xI_low = filter(b_low,a_low,xI_upsamp);
xQ_low = filter(b_low,a_low,xQ_upsamp);

% Compensate filter delay 
xI_low = xI_low((n_low/2)+1:end);
xQ_low = xQ_low((n_low/2)+1:end);

%----------------------Apply pulseform ------------------------------------

% Generate chirp-signal
t1 = 0:Ts:1-Ts;
f0 = 1;          %Start frequency
epsilon = 10;    %Chirp-rate
chirp = cos(2*pi*f0*(1+(epsilon.*t1)).*t1);

% Generate some extra zeros due to noise from the channel
extra_zeros = zeros(1, 40);
chirp = [chirp  extra_zeros];

% Apply chirp to xI_low and zeroes to xQ_low to match dimensions.
xI_with_chirp = [transpose(chirp); xI_low];
xQ_with_zeros = [zeros(length(chirp),1); xQ_low];

%------------------Generate signal to send via I/Q-------------------------

% I/Q - modulation (5 second signal, 1 second chirp)
t2 = 0:Ts:(length(xI_with_chirp))* Ts - Ts;
I_carry = transpose(cos(2*pi*fc*t2));
Q_carry = transpose(-sin(2*pi*fc*t2));
xI_send = xI_with_chirp .* I_carry;
xQ_send = xQ_with_zeros .* Q_carry;

x = xI_send + xQ_send;

end