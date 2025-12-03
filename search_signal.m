clear
fs = 1000;            % Частота дискретизации, Гц
duration = 100;         % Длительность сигнала, с
t = 0:1/fs:duration-1/fs; % временная ось
A = 3
signal1 = 10*sin(325*pi*2*t)
signal3 = 10 *sin(375*pi*2*t)
SNR = 59
f_start = 200;
f_end = 400;
signal2 = A*chirp(t, f_start, duration, f_end, 'linear');
f_start = 275;
f_end = 350;
signal4 = 15*chirp(t, f_start, duration, f_end, 'linear');
signal =  (sqrt(1/SNR)*A)*randn(size(t))+signal2+signal4+signal1+signal3;

nfft = 4096;
window = hamming(512);
noverlap = 256;

[spectrum, freq_axis] = pwelch(signal, window, noverlap, nfft, fs);
freq_range = [0 500];
freq_mask = freq_axis >= freq_range(1) & freq_axis <= freq_range(2);
freq_axis_limited = freq_axis(freq_mask);
spectrum_limited = spectrum(freq_mask);

window_size = 80;
threshold_db = 5;
detected_signals = detector_signals(freq_axis_limited, sqrt(spectrum_limited), window_size, threshold_db);

figure('Position', [100, 100, 1000, 800]);
subplot(2,1,1);
plot(freq_axis_limited, 10*log10(spectrum_limited), 'k', 'LineWidth', 1.2);
hold on; grid on;
xlabel('Частота (Гц)');
ylabel('Мощность (дБ)');
title('Спектр сигнала');
xlim(freq_range);
yl = ylim;

colors = lines(length(detected_signals));

for i = 1:length(detected_signals)
    sig = detected_signals(i);
    c = colors(i,:);
 
    line([sig.start_freq sig.start_freq], yl, 'Color', c*0.7 + 0.3, 'LineWidth', 1, 'LineStyle', '-');
    line([sig.end_freq sig.end_freq],     yl, 'Color', c*0.7 + 0.3, 'LineWidth', 1, 'LineStyle', '-');
    line([sig.center_freq sig.center_freq], yl, 'Color', c*0.7 + 0.3, 'LineWidth', 1, 'LineStyle', '--');
end

fprintf('Обнаруженные сигналы\n');
for i = 1:length(detected_signals)
    sig = detected_signals(i);
    fprintf('%d: %s, Центр=%.1f Гц, Ширина=%.1f Гц\n', i, ru_type(sig.type), sig.center_freq, sig.bandwidth);
end
% это для того, что бы в консоль вводился текст на русском языке
function str = ru_type(type_en)
    switch lower(type_en)
        case 'narrowband'
            str = 'Узкополосный';
        case 'broadband'
            str = 'Широкополосный';
        otherwise
            str = type_en;
    end
end
