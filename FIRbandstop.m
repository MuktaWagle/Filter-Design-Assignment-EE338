f_samp = 260e3;

%Band Edge speifications
fs1 = 55.1e3;
fp1 = 59.1e3;
fp2 = 79.1e3;
fs2 = 83.1e3;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-7.95) / (2.285*0.0242*pi));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 3;

%Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) - ideal_lp(0.6235*pi,n) + ideal_lp(0.439*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response

sl = [55.1e3 55.1e3 55.1e3];
pl = [59.1e3 59.1e3 59.1e3];
ph = [79.1e3 79.1e3 79.1e3];
sh = [83.1e3 83.1e3 83.1e3];
y = [0.0 0.75 1.2];
x = [0e4 10e4 18e4];
del1 = [1.15 1.15 1.15];
del2 = [0.85 0.85 0.85];
del3 = [0.15 0.15 0.15];

x = -26:26;
stem(x,FIR_BandStop,'.')
xlabel('n');
ylabel('h[n]');
xlim([-35,35]);
title('Time Domain Sequence');

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
plot(f,abs(H),'b',sl,y,pl,y,ph,y,sh,y,x,del1,x,del2,x,del3)
xlim([0, 13e4]);
grid
