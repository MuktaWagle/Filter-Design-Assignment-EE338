f_samp = 330e3;

%Band Edge speifications
fs1 = 66.1e3;
fp1 = 70.1e3;
fp2 = 90.1e3;
fs2 = 94.1e3;

Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-7.95) / (2.285*0.0242*pi));           %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 13;

%Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(0.558*pi,n) - ideal_lp(0.4125*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass);         %frequency response

sl = [66.1e3 66.1e3 66.1e3];
pl = [70.1e3 70.1e3 70.1e3];
ph = [90.1e3 90.1e3 90.1e3];
sh = [94.1e3 94.1e3 94.1e3];
y = [0.0 0.75 1.2];
x = [0e4 10e4 18e4];
del1 = [1.15 1.15 1.15];
del2 = [0.85 0.85 0.85];
del3 = [0.15 0.15 0.15];

x1 = -31:31;
stem(x1,FIR_BandPass,'.') %time domain sequence plot
xlabel('n');
ylabel('h[n]');
xlim([-35,35]);
title('Time Domain Sequence'); 

%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, f_samp);
plot(f,abs(H),'b',sl,y,pl,y,ph,y,sh,y,x,del1,x,del2,x,del3);
grid on
grid minor
