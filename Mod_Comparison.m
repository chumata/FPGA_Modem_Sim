clear variables
close all
clc

cd('D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim')

%Open txt file
file_name='D:\Projects\SDRZ\Simulations\FPGA_Modem_Sim\VHDL_IQ_Output\qam\16\QAM16.txt';
data=open_txt(file_name);

M=16; %symbol order
sam_num=50; %Number of samples per symbol

[data_VHDL,QAM_sim,QAM_diff]=QAM_Comparison(M,sam_num,data);

function data=open_txt(file_name)

    fid =fopen(file_name);
    data = textscan(fid,'%f%f','HeaderLines',2,'CollectOutput',1);
    data = data{:};
    fclose(fid);
end

function [data_VHDL,QAM_sim,QAM_diff]=QAM_Comparison(M,sam_num,data)
    x=(0:M-1)'; %Symbols vector
    I_Q= qammod(x,M); %Grey code QAM IQ
    t=0:1/sam_num:1-1/sam_num; %time vector
    QAM_sim=zeros(sam_num,M);
    for i=1:M
        QAM_sim(:,i)=real(I_Q(i))*cos(2*pi*t+pi/2)-imag(I_Q(i))*sin(2*pi*t+pi/2);
    end
    QAM_sim=reshape(QAM_sim,[],1);
    power_coeff=max(data(:,2))/max(QAM_sim); %Normalize power
    data_VHDL=data(:,2)/power_coeff;
    QAM_sim=QAM_sim(1:length(data_VHDL));
    QAM_diff=QAM_sim-data_VHDL;
    figure
    subplot(2,1,1)
    plot(data_VHDL)
    hold on
    plot(QAM_sim)
    grid on
    legend('VHDL','MATLAB')
    title('VHDL vs. matlab comparison for XXX modulation')
    xlabel('samples')
    ylabel('amplitude')
    subplot(2,1,2)
    plot(QAM_diff)
    grid on
    title('VHDL vs. matlab difference')
    xlabel('samples')
    ylabel('amplitude')

end