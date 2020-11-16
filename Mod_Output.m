clear variables
close all
clc
cd('D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim')
%%
%Symbol generator for FPGA correlator.
%Cases:
%mod=1 - BPSK- declare modulation order M=2.
%mod=2 - QPSK- declare modulation order M=4.
%mod=3 - QAM- declare modulation order with variable M. Different branches I-Q.
%mod=4 - QAM- declare modulation order with variable M. Single branch I+Q.


mod=3; %case
bit_num=13; %Number of bits representing signed amplitude
N=10; %Number of samples per symbol

switch mod
    case 1
    %%BPSK
    M=2; %modulation order
    path = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\BPSK\';  
    sym=int16(BPSK_symbols(bit_num,N,M)); %Integer rows vector modulated symbols
    write_txt(M,sym,path); %Write and save formatted text files
    
    case 2
    %%QPSK
    M=4; %modulation order
    path_I = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QPSK\I\';
    path_Q = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QPSK\Q\';  
    [I_sym,Q_sym]=QPSK_symbols(bit_num,N,M); %Rows vector modulated symbols
    int_I_sym=int16(I_sym); %Integer rounding
    int_Q_sym=int16(Q_sym); %Integer rounding
    write_txt(M,int_I_sym,path_I); %Write and save formatted text files
    write_txt(M,int_Q_sym,path_Q); %Write and save formatted text files
    
    case 3
    %%QAM
    M=16; %modulation order
    [I_sym,Q_sym]=QAM_symbols(bit_num,N,M);
    int_I_sym=int16(I_sym); %Integer rounding
    en_int_I_sym=sum(int32(I_sym).^2,2); %symbol energy
    I_TH=en_int_I_sym(sqrt(M)+1)+(en_int_I_sym(sqrt(M))-en_int_I_sym(sqrt(M)+1))/2;
    int_Q_sym=int16(Q_sym); %Integer rounding
    en_int_Q_sym=sum(int32(Q_sym).^2,2); %symbol energy
    Q_TH=en_int_Q_sym(2)+(en_int_Q_sym(1)-en_int_Q_sym(2))/2;
    switch M
        case 4
        path_I = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QAM\4QAM\I\';
        path_Q = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QAM\4QAM\Q\'; 
        write_txt(M,int_I_sym,path_I); %Write and save formatted text files
        write_txt(M,int_Q_sym,path_Q); %Write and save formatted text files
        case 16
        path_I = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QAM\16QAM\I\';
        path_Q = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QAM\16QAM\Q\'; 
        write_txt(M,int_I_sym,path_I); %Write and save formatted text files
        write_txt(M,int_Q_sym,path_Q); %Write and save formatted text files
        case 64
        path_I = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QAM\64QAM\I\';
        path_Q = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QAM\64QAM\Q\'; 
        write_txt(M,int_I_sym,path_I); %Write and save formatted text files
        write_txt(M,int_Q_sym,path_Q); %Write and save formatted text files
    end
    
        case 4
    %%QAM IQ
    M=16; %modulation order
    IQ_sym=QAM_symbols_I_minus_Q(bit_num,N,M);
    int_IQ_sym=int16(IQ_sym); %Integer rounding
    en_int_IQ_sym=sum(int32(IQ_sym).^2,2); %symbol energy
    switch M
        case 4
        path_IQ = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QAM\4QAM\I_minus_Q\'; 
        write_txt(M,int_IQ_sym,path_IQ); %Write and save formatted text files
        case 16
        path_IQ = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QAM\16QAM\I_minus_Q\'; 
        write_txt(M,int_IQ_sym,path_IQ); %Write and save formatted text files
        case 64
        path_IQ = 'D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\Symbols_Output_Sim\QAM\64QAM\I_minus_Q\';
        write_txt(M,int_IQ_sym,path_IQ); %Write and save formatted text files
    end
end

function sym=BPSK_symbols(bit_num,N,M)
    t=linspace(0,1,N+1); %time vector
    A=2^bit_num-1; %Amplitude (normalized for 14 bit signed fixed point)
    sym_vec=0:1:M-1; %symbols vector
    bpskModulator = comm.BPSKModulator;
    modData = bpskModulator(sym_vec');
    sig=A*real(modData).*cos(2*pi*t);
    sym=sig(:,1:end-1);
    figure
    scatterplot(modData)
    grid on
    figure
    plot(sym(1,:))
    grid on
    hold on
    plot(sym(2,:))
    title('BPSK')
end

function [I_sym,Q_sym]=QPSK_symbols(bit_num,N,M)
    t=linspace(0,1,N+1); %time vector
    A=2^bit_num-1; %Amplitude (normalized for 14 bit signed fixed point)
    sym_vec=0:1:M-1; %symbols vector
    qpskModulator = comm.QPSKModulator;
    modData = qpskModulator(sym_vec');
    I_sig=A*(real(modData)).*cos(2*pi*t);
    I_sym=I_sig(:,1:end-1);
    Q_sig=A*(imag(modData)).*sin(2*pi*t);
    Q_sym=Q_sig(:,1:end-1);
    figure
    scatterplot(modData)
    grid on
    figure
    %Plot I channel symbols
    for i=1:M
    plot(I_sym(i,:))
    grid on
    hold on
    title('I Channel QPSK')
    end
    figure
    %Plot Q channel symbols
    for i=1:M
    plot(Q_sym(i,:))
    grid on
    hold on
    title('Q Channel QPSK')
    end   
end

function [I_sym,Q_sym]=QAM_symbols(bit_num,N,M)
    t=linspace(0,1,N+1); %time vector
    A=2^bit_num-1; %Amplitude (normalized for 14 bit signed fixed point)
    sym_vec=0:1:M-1; %symbols vector
    modData = (qammod(sym_vec,M,'UnitAveragePower',true,'PlotConstellation',true))';
    I_sig=A*(real(modData)).*cos(2*pi*t);
    I_sym=I_sig(:,1:end-1);
    Q_sig=A*(imag(modData)).*sin(2*pi*t);
    Q_sym=Q_sig(:,1:end-1);
    figure
    %Plot I channel symbols
    for i=1:M
    plot(I_sym(i,:))
    grid on
    hold on
    title('I Channel QAM')
    end
    figure
    %Plot Q channel symbols
    for i=1:M
    plot(Q_sym(i,:))
    grid on
    hold on
    title('Q Channel QAM')
    end   
end

function IQ_sym=QAM_symbols_I_minus_Q(bit_num,N,M)
    t=linspace(0,1,N+1); %time vector
    A=2^bit_num-1; %Amplitude (normalized for 14 bit signed fixed point)
    sym_vec=0:1:M-1; %symbols vector
    modData = (qammod(sym_vec,M,'UnitAveragePower',true,'PlotConstellation',true))';
%     I_sig=A*(real(modData)).*cos(2*pi*t);
%     I_sym=I_sig(:,1:end-1);
%     Q_sig=A*(imag(modData)).*sin(2*pi*t);
%     Q_sym=Q_sig(:,1:end-1);
    IQ_sig=A*((real(modData)).*cos(2*pi*t)-(imag(modData)).*sin(2*pi*t));
    IQ_sym=IQ_sig(:,1:end-1);
    figure
    %Plot I channel symbols
    for i=1:M
    plot(IQ_sym(i,:))
    grid on
    hold on
    title('I-Q Channel QAM')
    end
end

function write_txt(M,sig,path)
pre="to_signed(";
after=",14),";
vec=strings(M,1);
for m=1:M
str = string(sig(m,:));
    for i=1:length(str)
    newStr = strcat(pre,str(i),after);
    vec(m)=strcat(vec(m),newStr);
    end
end

for i=1:length(vec)
    T = table(vec(i));
    sym_name = sprintf("%d.txt",i-1);
    
    file_name = strcat(path,sym_name);
    writetable(T,file_name,'WriteVariableNames',0)
end
end

