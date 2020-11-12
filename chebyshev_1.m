%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 4;

% Open CLHP Poles of the Chebyshev Polynomial of order 4
B = asinh(1/epsilon)/N ;
p1 = -sin(pi/(2*N))*sinh(B)+1i*cos(pi/(2*N))*cosh(B);
p2 = -sin(pi/(2*N))*sinh(B)-1i*cos(pi/(2*N))*cosh(B);
p3 = -sin(3*pi/(2*N))*sinh(B)+1i*cos(3*pi/(2*N))*cosh(B);
p4 = -sin(3*pi/(2*N))*sinh(B)-1i*cos(3*pi/(2*N))*cosh(B);        

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2); %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(5)*sqrt(1/(1+D1))];        % even order, DC Gain set as 1/(1+D1)^0.5

%Band Edge speification
fs1 = 66.1;
fp1 = 70.1;
fp2 = 90.1;
fs2 = 94.1;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 330;
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bpf(s) = analog_lpf((s*s +W0*W0)/(B*s));     %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BPF
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete BPF
[nz, dz] = numden(discrete_bpf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                       %frequency response in dB

sl = [66.1e3 66.1e3 66.1e3];
pl = [70.1e3 70.1e3 70.1e3];
ph = [90.1e3 90.1e3 90.1e3];
sh = [94.1e3 94.1e3 94.1e3];
y = [0.0 0.75 1.2];
x = [0e4 10e4 18e4];
del1 = [1.15 1.15 1.15];
del2 = [0.85 0.85 0.85];
del3 = [0.15 0.15 0.15];

zplane(nz,dz);

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024,330e3);
plot(f,abs(H),'b',sl,y,pl,y,ph,y,sh,y,x,del1,x,del2,x,del3)
grid on 
grid minor