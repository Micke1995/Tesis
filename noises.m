clear all
close all
clc
%% Signals
Fs = 20;
Ts = 1/Fs;
Fmin = 0.2;
T1 = 1/Fmin;
N = round(T1*Fs);
n = 0:2*N-2;
t = n;
F1 = 0.23;
F2 = 0.36;
F3 = 0.64;
sigma1 = 0.01;
sigma2 = 0.001;
sigma3 = 0.05;
A1 = 2;
A2 = 3;
A3 = 1;
ph1 = pi/2;
ph2 = pi/3;
ph3 = pi/5;
%
s1 = A1*exp(-sigma1*n).*sin(2*pi*F1*n/Fs + ph1);
s2 = A2*exp(-sigma2*n).*sin(2*pi*F2*n/Fs + ph2);
s3 = A3*exp(-sigma3*n).*sin(2*pi*F3*n/Fs + ph3);
%
signals = s1 + s2 + s3;
[ns, Nm] = size(signals);
st = zeros(ns, Nm);
means = zeros(ns, Nm);
ra = zeros(ns, Nm);
snra = zeros(ns, 1);
%% 15 db noise signal
p = 23;%[96 97 97 92 92 92 93 96 93 99 100 98 97 105 103 104];
for k = 1:ns
    means(k,:) = mean(signals(k,:));
    st(k,:) = signals(k,:) - means(k,:);
    ra(k,:) = wgn(Nm, 1, -p(k));
    snra(k,1) = 10*log10( sqrt(mean(st(k,:).^2)) / sqrt(mean(ra(k,:).^2)));
    mean(st(k,:).^2)
end
stna = st + ra + means;
%
figure(1)
plot(t, signals, t, stna)
% save 15dB_angles_16m_vmd t mac_spd_15dB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%El programa anterior tiene los siguientes errores:
%a). La energ??a del ruido (-10dB) es independiente de la energ??a de
%la se??al, cuando deber??a de ser dependiente
%b). La funci??n wgn entrega un ruido con VARIANCIA de -10dB cuando queremos
%un ruido de ENERGIA. Se tiene En=(N-1)*sigma^2, N el n??mero de muestras
%b). El SNR en la variable snra deber??a de ser 20*log10() y no 10*log10()
%porque el argumento lleva ra??z cuadrada.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% El programa deber??a de ser como sigue para 20dB de SNR:
% SNRdB = 100;
% N = Nm;
% for k = 1:ns
%     Es(k,1) = sum(signals(k,:).^2);
%     EsdB(k,1) = 10*log10(Es(k,1));
%     EndB(k,1) = EsdB(k,1) - SNRdB;
%     vardB(k,1) = EndB(k,1) - 10*log10(N-1);
%     r(k,:) = wgn(Nm, 1, vardB(k,1))';
%     snrdB(k,1) = 10*log10( sum(signals(k,1).^2) / sum(r(k,:).^2) );  % SNR en decibeles
% end
% mac_spd_66dB = signals + r;
% 
% figure(2)
% plot(t, signals, t, mac_spd_66dB)
% 
% noise = [snra snrdB];
% disp(noise)
% % save job_68_fault_32n_46dB t mac_spd_46dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Como ven en la figura 2 el ruido de SNR=20dB es mucho m??s grande que el
%que ellos afirman que tiene 21.0346
%Comparando las graficas, el ruido de ellos es equivalente a un ruido de
%40dBs (corran mi programa con SNRdB=40 y dan visualmente parecidos.
%Si se centrara la se??al, la varianza ser??a proporcional a la energ??a y la
%variancia del ruido tambi??n y los factores 1/(N-1) de ambas variancias se
%cancelar??an en sigma_S^2/sigma_N^2=Es/En, pero como la energ??a de la se??al
%es sin centrar, hay que proceder como el parche que les env??o.
