%Butterworth Analog LPF parameters
Wc = 1.075;              %cut-off frequency
N = 7;                  %order 

%poles of Butterworth polynomial of degree 7 in the open CLHP 
p1 = Wc*cos(pi/2 + pi/14) + 1i*Wc*sin(pi/2 + pi/14);
p2 = Wc*cos(pi/2 + pi/14) - 1i*Wc*sin(pi/2 + pi/14);
p3 = Wc*cos(pi/2 + pi/14+ pi/7) + 1i*Wc*sin(pi/2 + pi/14+ pi/7);
p4 = Wc*cos(pi/2 + pi/14+ pi/7) - 1i*Wc*sin(pi/2 + pi/14+ pi/7);
p5 = Wc*cos(pi/2 + pi/14+ 2*pi/7) + 1i*Wc*sin(pi/2 + pi/14+ 2*pi/7);
p6 = Wc*cos(pi/2 + pi/14+ 2*pi/7) - 1i*Wc*sin(pi/2 + pi/14+ 2*pi/7);
p7 = Wc*cos(pi/2 + pi/14+ 3*pi/7) + 1i*Wc*sin(pi/2 + pi/14+ 3*pi/7);


%Band Edge speifications
fp1 = 55.1;
fs1 = 59.1;
fs2 = 79.1;
fp2 = 83.1;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 260;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                           %frequency response

sl = [55.1e3 55.1e3 55.1e3];
pl = [59.1e3 59.1e3 59.1e3];
ph = [79.1e3 79.1e3 79.1e3];
sh = [83.1e3 83.1e3 83.1e3];
y = [0.0 0.75 1.2];
x = [0e4 10e4 12.99e4];
del1 = [1.15 1.15 1.15];
del2 = [0.85 0.85 0.85];
del3 = [0.15 0.15 0.15];

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 260e3);
plot(f,abs(H),'b',sl,y,pl,y,ph,y,sh,y,x,del1,x,del2,x,del3)
grid