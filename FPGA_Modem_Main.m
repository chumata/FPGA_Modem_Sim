clear variables
close all
clc

cd('D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim')

A=1; %Sine wave amplitude
Fc=33.25e6; %Carrier frequency [Hz]
T_sym=1e-6; %Symbol period time [seconds]
Sym_rate=1/T_sym; %Symbol rate [sym/sec]
Mod_Order=4; %Modulation order [2,4,8...]
T_bit=T_sym/log2(Mod_Order); %Bit period time [seconds]
Bit_rate=1/T_bit; %Bit rate [bps]
Fs=150e6; %Continious wave frequency
Bits_Num=8; %A2D Number of bits
A2D_Fs=60e6; %A2D sampling frequency [samples/seconds]
Frame_Size=64; %Time samples to be processed. The maximum for FFT block input is 64.
alpha = 0.8; %A2D non- linearity factor
FFT_Length=512; %FFT resolution 
RBW=Fs/FFT_Length; %Frequency bin=Spectrum RBW
Span=FFT_Length*RBW/2; %Frequency span (the FFT output is double sided)
Quantizer_Mapping=0:2^-Bits_Num:1; %Quantizer options
Frame_Periods=1000;
stop_time=Frame_Periods*Frame_Size/Fs; %Simulation stop time
jitter=5*1e-12; %Max clock jitter [seconds]
Freq_Stable=0.5*1e-6; %Freq
Dec_factor_1=30; %1st stage decimation factor
Phase_Offset=0; %Signal phase offset
fix_WL=10; %Fixed point word length
fix_FL=8; %Fixed point fraction length out of the word length
beta_RC=0.2; %Square root cosine filter roll off
L_RC=4; %Pulse shaping filter interpolation coefficient
%% open Simulink model
sim_file='FPGA_Modem_Sim.slx';
open_system(sim_file);