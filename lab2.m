clc
clear 
close all

%% формирование сетки

Nfft = 2048;
Gi = 1/16;

all_sc = -Nfft/2 : Nfft/2 - 1; 
all_sc_shifted = all_sc+Nfft/2+1;
zero_sc = Nfft/2+1;

Gi_left = 1:Nfft*Gi/2+1;
Gi_right = (Nfft-Nfft*Gi/2+1):Nfft;

pilot_sc = max(Gi_left):16:min(Gi_right)-1;

exc_sc = [Gi_left,Gi_right,zero_sc,pilot_sc];%pilot_sc
Y_inf = all_sc_shifted(~ismember(all_sc_shifted,exc_sc));
%% Формирование парметров битов
M_inf = 16;
k_inf = log2(M_inf);
n_inf_bits = length(Y_inf);
M_pilot = 4;
k_pilot = log2(M_pilot);
n_pilot_bits = length(pilot_sc);
%% формирование пилотов
rng(1);
pilot_bits = randi([0,1],n_pilot_bits*k_pilot,10);
qampilot = qammod(pilot_bits,M_pilot,'InputType', 'bit','UnitAveragePower', true);

%% Константы кадра
CP_len = Nfft/4;%Nfft/4;
T_symb = 100e-6;
Fs = Nfft/T_symb;
%% Канал EPA/EVA/ETU

path_delay_EPA = [0 30 70 80 110 190 410]*10^-9;
APG_EPA = [0 -1 -2 -3 -8 -17.2 -20.8];
rayleigh_EPA = comm.RayleighChannel('SampleRate',Fs,'PathDelays',path_delay_EPA,'AveragePathGains',APG_EPA, 'MaximumDopplerShift', 0,'NormalizePathGains',true);

path_delay_EVA = [0 30 150 310 370 710 1090 1730 2510]*10^-9;
APG_EVA = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9];
rayleigh_EVA = comm.RayleighChannel('SampleRate',Fs,'PathDelays',path_delay_EVA,'AveragePathGains',APG_EVA, 'MaximumDopplerShift', 0,'NormalizePathGains',true);

path_delay_ETU = [0 50 120 200 230 500 1600 2300 5000]*10^-9;
APG_ETU = [-1 -1 -1 0 0 0 -3 -5 -7];
rayleigh_ETU = comm.RayleighChannel('SampleRate',Fs,'PathDelays',path_delay_ETU,'AveragePathGains',APG_ETU, 'MaximumDopplerShift', 0,'NormalizePathGains',true);
%% ИХ и ЧХ каналов

Nfft_ch = 2^12;
Fs_ch = 10^9;
T_ch = 1/Fs_ch;
df = Fs_ch/Nfft_ch;
f = -Fs_ch/2:df:(Fs_ch/2-df);

tao = path_delay_EPA;
P_db = APG_EPA;
P = 10.^(P_db / 20);

h = zeros(1,round(max(tao)*Fs_ch+1));
h(round(tao*Fs_ch+1)) = sqrt(P);

H = fftshift(fft(h,Nfft_ch));
%графики
figure 
plot((0:length(h)-1)*T_ch, abs(h));
figure 
plot(f,abs(H));
xlim([-12e6, 12e6]);
%% канал 
SNR_vector = 0:30; 
BER_vector = zeros(1, length(SNR_vector));
N_frames = 100;   

for idx_ber = 1:length(SNR_vector)
    SNR = SNR_vector(idx_ber);
%     SNR = 30;
    
    total_errors = 0;
    total_bits = 0;
    
    for frame = 1:100 
        
        bits = randi([0,1],n_inf_bits*k_inf,10);

        qamsymb = qammod(bits,M_inf,'InputType', 'bit','UnitAveragePower', true);
        % Формируем матрицу переданных бит в нужном формате
        bits_dem = qamdemod(qamsymb.', 16,'OutputType', 'bit','UnitAveragePower', true);

        % формирование OFDM символа
        tx_fd_sig = zeros(10,Nfft);
        tx_fd_sig(:,pilot_sc) = qampilot.';
        tx_fd_sig(:,Y_inf) = qamsymb.';
        OFDM_symb = ifft(ifftshift(tx_fd_sig,2),[],2);

        % добавление ЦП
        OFDM_symb_cp = [OFDM_symb(:,end-CP_len+1:end),OFDM_symb];%zeros(10,CP_len) или OFDM_symb(:,end-CP_len+1:end)

        % АБГШ       
        SNR_linear = 10^(SNR / 10);
        
        % Измеряем точную среднюю мощность чистого сигнала до добавления ЦП
        P_sig = mean(abs(OFDM_symb).^2,'all'); 
        
        % Рассчитываем необходимую мощность шума
        P_noise = P_sig / SNR_linear;
        
        % Амплитуда комплексного шума
        sigma = sqrt(P_noise / 2);
        
        % Генерируем комплексный гауссовский шум
        noise = sigma * (randn(size(OFDM_symb_cp)) + 1i * randn(size(OFDM_symb_cp)));
        %преобразование в один сигнал для передачи
        OFDM_vector = OFDM_symb_cp.';
        OFDM_vector = OFDM_vector(:);

        % релей
        OFDM_rel = rayleigh_EPA(OFDM_vector);
    
        OFDM_noise = reshape(OFDM_rel,size(OFDM_symb_cp,2),size(OFDM_symb_cp,1));% если без релея то OFDM_vector, иначе OFDM_rel
        OFDM_noise = OFDM_noise.'+noise;
        % приемник

        % удаление цп
        OFDM_rx = OFDM_noise(:,CP_len+1:end);
        
        % переход в частоты
        OFDM_rx_fd = fftshift(fft(OFDM_rx,[],2),2);
        
        % Извлечение данных
        rx_symbs = OFDM_rx_fd(:,Y_inf);
%         scatterplot(rx_symbs(1,:));
        rx_pilots = OFDM_rx_fd(:,pilot_sc);
        H = rx_pilots./tx_fd_sig(:,pilot_sc);
        
        H_all = interp1(pilot_sc,mean(H),1:Nfft,"linear","extrap");

        rx_data_equalized = OFDM_rx_fd./H_all;
        rx_data_equalized = rx_data_equalized(:,Y_inf); 
        rx_data_demap = qamdemod(rx_data_equalized, 16,'OutputType', 'bit','UnitAveragePower', true);% если без релея то = rx_symbs, с релеем rx_data_equalized
        % Считаем ошибки текущего фрейма
        [numErrors, ~] = biterr(bits_dem(:), rx_data_demap(:));

        total_errors = total_errors + numErrors;
        total_bits = total_bits + numel(bits);
%         scatterplot(rx_data_equalized(1,:));
        
    end
%     Сохраняем усредненный BER для текущего SNR
    BER_vector(idx_ber) = total_errors / total_bits;
end



%% Теоретическая кривая

% Коэффициент частотной избыточности
L_data = length(Y_inf) / Nfft; 

EbN0_vector = SNR_vector-10*log10(k_inf)-10*log10(L_data);
EbN0_linear = 10.^(EbN0_vector/10);

% Теоретический BER для 16-QAM в канале АБГШ
BER_theory = berawgn(EbN0_vector, 'qam', M_inf);

BER_theory_Rayleigh = berfading(EbN0_vector, 'qam', M_inf,1);
%% Построение графика

figure;

% 1. текущий эксперимент в канале с замираниями (Синий маркер)
semilogy(SNR_vector, BER_vector, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;

% 2. Теория для канала замираний Рэлея (Зеленый пунктир) — ЭТАЛОН ДЛЯ ЭТОГО ТЕСТА
semilogy(SNR_vector, BER_theory_Rayleigh, 'g--', 'LineWidth', 2);

% 3. Теория для канала АБГШ (Красный точечный) — ваша проверенная идеальная прямая
semilogy(SNR_vector, BER_theory, 'r:', 'LineWidth', 1.5);

hold off;
grid on;

% Настройка масштаба осей
xlim([min(SNR_vector), max(SNR_vector)]);
ylim([1e-4, 1]); 

xlabel('SNR, дБ');
ylabel('BER');
title('Кривая помехоустойчивости OFDM-системы в канале с замираниями');

legend('Экспериментальный BER (с эквалайзером)', ...
       'Теоретический BER (Рэлей, Эталон)', ...
       'Теоретический BER (Идеальный АБГШ)', ...
       'Location', 'SouthWest');
%% График для сравнения бер каналов
figure;


semilogy(SNR_vector, EPA_BER, 'b-o', 'LineWidth', 1, 'MarkerFaceColor', 'b');
hold on;


semilogy(SNR_vector, EVA_BER, 'g--', 'LineWidth', 1);


semilogy(SNR_vector, ETU_BER, 'r:', 'LineWidth', 1);

hold off;
grid on;

% Настройка масштаба осей
xlim([min(SNR_vector), max(SNR_vector)]);


xlabel('SNR, дБ');
ylabel('BER');
title('Кривая помехоустойчивости OFDM-системы в канале с замираниями');

legend('EPA BER', ...
       'EVA BER', ...
       'ETU BER', ...
       'Location', 'SouthWest');