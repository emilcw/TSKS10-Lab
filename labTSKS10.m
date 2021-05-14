% Laboration i kursen TSKS10
% Lösning av Emil Wiman, VT 2021
% Här finns main-scriptet som kallar på
% hjälpfunktioner i andra filer för att lösa
% uppgiften. 


%Tilldelning av audiosignaler
[xI,fs] = audioread('xI.wav');
[xQ,fs] = audioread('xQ.wav');

% Ta fram insignal till kanal, x kolumnvektor
x = sender(xI,xQ);

%Skicka via kanal och få ut utsignal
y = TSKS10channel(x);

%Hämta ut ny signal och skalning samt fördröjning
[zI,zQ,A,tau] = receiver(y);
A
tau
%figure(1)
%plot(xI)
%title('xI');
%figure(2)
%plot(zI)
%title('zI')
%figure(3)
%plot(xQ)
%title('xQ')
%figure(4)
%plot(zQ)
%title('zQ')

% Spela upp ny signal
%soundsc(zI, fs)
%soundsc(zQ, fs);

% Kontroll av utsignal, svar bör vara minst 25 db
SNRzI = 20*log10(norm(xI)/norm(zI-xI))
SNRzQ = 20*log10(norm(xQ)/norm(zQ-xQ))

