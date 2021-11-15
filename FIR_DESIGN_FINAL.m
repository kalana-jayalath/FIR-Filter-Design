clear all;
close all ;
%M.W.K.S.Jayalath
%Index number - 180260U

%taking the values needed for A , B , C extracted from the index number
Index_number = input('Enter your index number:');
num_A = mod(floor(Index_number/100),10);
num_B = mod(floor(Index_number/10),10);
num_C = mod(floor(Index_number),10);

%disp(class(num_C))
%disp(num_C)

%Calculation of parameters
Max_Ap = 0.03 + (0.01*num_A);            %Maximum passband ripple, A˜p 0.03 + (0.01 × A) dB
Min_Aa = (45 + num_B);                   %Minimum stopband attenuation, A˜a 45+B dB
Omega_LP = (num_C*100)+300 ;              %Lower passband edge, 
                                        %p1 (C × 100) + 300 rad/s
Omega_UP = (num_C*100)+700  ;            %Upper passband edge, 
                                        %p2 (C × 100) + 700 rad/s
Omega_LS = (num_C*100)+150  ;            %Lower stopband edge, 
                                        %a1 (C × 100) + 150 rad/s
Omega_US = (num_C*100)+800  ;            %Upper stopband edge, 
                                        %a2 (C × 100) + 800 rad/s
Omega_smp = 2*((num_C*100)+1200);        %Sampling frequency, 
                                        %s 2[(C × 100) + 1200] rad/s

%calculation of the delta value and stop band attenuation
del_p = (10^(0.05*Max_Ap)-1) / (10^(0.05*Max_Ap)+1) ; 
del_a = (10^(-0.05*Min_Aa));
del   = min(del_p , del_a);
stop_at = -20*log10(del);

%calculation of cut off frequencies and transition bandwidth
B_t = min( (Omega_US - Omega_UP) , (Omega_LP - Omega_LS) );
cf_up = (Omega_UP+(B_t/2));
cf_lw = (Omega_LP-(B_t/2));

%calculation of alpha
if stop_at < 21
    alpha = 0;
elseif 21<stop_at && stop_at <= 50
    alpha = 0.5842*(stop_at-21)^0.4 + 0.07886*(stop_at-21) ;
else
    alpha = 0.1102*(stop_at-8.7) ;
end

%calculation of D
if stop_at <= 21
    D = 0.9222;
else 
    D = (stop_at - 7.95)/ 14.36;
end

%calculation of N
if mod( ceil(Omega_smp*D/B_t + 1) ,2) == 1
    N = ceil(Omega_smp*D/B_t + 1) ;
else
    N = ceil(Omega_smp*D/B_t + 1) + 1;
end

disp(N)

%computation and plotting the kaizer window
half_range = (N-1) /2;
n = -half_range : 1 : half_range ;
beta =alpha * (1 - ((2*n/(N-1)).^2)).^0.5;
I_beta = 0; 
I_alpha = 0;
for k = 1 : 100
    I_beta = I_beta + ((1/factorial(k))*(beta/2).^k).^2;
    I_alpha = I_alpha + ((1/factorial(k))*(alpha/2)^k)^2;
end

I_beta = I_beta + ones(1 ,numel( I_beta ) ) ;
I_alpha = I_alpha + ones(1 ,numel( I_alpha ) ) ;
w = I_beta/I_alpha ;

figure ; 
stem(n,w,'fill' ) ;
xlabel('n') ;
ylabel ('w[n]') ; 
title( 'kaizer window') ;
grid on;


%computing ideal bandpass filter impulse respose
T=2*pi/Omega_smp;
half_range=(N-1)/2;

n_neg=-half_range:1:-1;
hI_neg=  ((1./n_neg)*(1/pi)).*(sin(cf_up*n_neg*T)-sin(cf_lw*n_neg*T));

hI_0=2*(cf_up-cf_lw)/Omega_smp;

n_pos=1:1:half_range;
hI_pos= ((1./n_pos)*(1/pi)).*(sin(cf_up*n_pos*T)-sin(cf_lw*n_pos*T));


hI=[hI_neg,hI_0,hI_pos];
n=[n_neg,0,n_pos];

figure;
stem(n,hI,'fill');
xlabel('n');
ylabel('hI[n]');
title('Ideal Impulse Response of a Bandpass Filter');
grid on;


%Computing impulse response for FIR bandpass filter
h_res=hI.*w; 

figure;
stem(n,hI,'fill');
xlabel('n');
ylabel('h[n]');
title('Impulse Response of a FIR Causal Bandpass Filter');
grid on;

shift_n=[0:1:N-1];
figure;
stem(shift_n,h_res,'fill');
xlabel('n');
ylabel('hI[n]');
title('Impulse Response of a FIR non Causal Bandpass Filte180r');
grid on;

%Magnitude response of the Bandpass Filter
fvtool(h_res);
freqz(h_res);

%Output of filter to a sudden excitation
Omega_one=Omega_LS/2;
Omega_two=(Omega_LP+Omega_UP)/2;
Omega_three=(Omega_US+Omega_smp)/2;

n=0:1:300;
m=length(n);
x=sin(Omega_one*n*T)+sin(Omega_two*n*T)+sin(Omega_three*n*T);
L=length(h_res);

%plotting DFT of the excitation signal
L = numel(n);
n_point = 2^nextpow2(L) ; % Next power of 2 from length of y
Y = fft(x ,n_point)/L;
f = (Omega_smp) /2* linspace (0 ,1 ,n_point/2+1) ;
figure ;
subplot (3 ,1 ,1) ;
plot( f ,2*abs (Y( 1 :n_point/2+1) ) ) ;
title('DFT of the excitation signal')
xlabel ('frequency in radians per second');
ylabel ( '|x(f)|' ) ; 
grid on;

%Plot the DFT of filtered signal

x_f = conv ( x , h_res ,'same' ) ;
L = numel ( x_f ) ;
n_point = 2^nextpow2(L) ; % Next power of 2 from length of y
Y = fft( x_f ,n_point) /L ;
f = (Omega_smp) /2* linspace (0 ,1 ,n_point/2+1) ;
subplot (3 ,1 ,2) ;
plot( f ,2* abs (Y( 1 :n_point/2+1) ) ) ;
title('DFT of the filtered signal' );
xlabel ('Frequency ( rad/ s )' );
ylabel (' |X( f ) | ') ; 
grid on;

%Plot the DFT of the excitation passed through an ideal bandpass fiter

x_i =conv ( x , hI ,'same' );
L = numel (n) ;
n_point = 2^nextpow2(L) ; % Next power of 2 from length of y
Y = fft( x_i ,n_point) /L ;
f = (Omega_smp) /2* linspace (0 ,1 ,n_point/2+1) ;
subplot (3 ,1 ,3) ;
plot( f ,2* abs (Y( 1 :n_point/2+1) ) ) ;
title('DFT of excitation that passed through an ideal bandpass filter') ;
xlabel ('Frequency ( rad/ s )')
ylabel (' |X( f ) | ') ;
grid on;

%Plot the excitation in time domain

figure ;
stem(n, x ,'r','fill' ) ;
title('Excitation in time domain')
xlabel ('n')
ylabel ('x[n]') ; 
grid on;

%Plot the Fi l tered s ignal in the time domain

nl = [0 : 1 :numel( x_f )-1];

figure ;
stem( nl , x_f , 'r', 'fill' ) ;
title( 'filtered sgnal in time domain');
xlabel ('n' )
ylabel ('x[n]') ; 
grid on;

%Plot the Ideal ly Fi l tered s ignal in the time domain

figure ;
stem(n, x_i ,'r','fill' ) ;
title('signal filtered using the ideal filter in time domain');
xlabel ('n');
ylabel ('x[n]') ; 
grid on;
    















