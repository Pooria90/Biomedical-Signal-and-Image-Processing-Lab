%% ======================================================= Part 1: Display
data    = load('data/EEG_sig.mat');
spec    = data.des;
sig     = data.Z;
fs      = spec.samplingfreq;

f1 = figure;
f1.WindowState = 'maximized';
%f1.Position = [1 41 1280 606];
% figure('units','normalized','outerposition',[0 0 1 1])
ch  = sig(5,:);
T   = (0:length(ch)-1)/fs;
plot(T,ch)
grid on
xlabel('Time (s)')
ylabel('Amplitude')
title(['Signal on channel ',spec.channelnames{5}])


%% ==================================================== Part 2: Partitions
I_1 = find(T<15);
T_1 = T(I_1);
P_1 = ch(I_1);

I_2 = find(T>18 & T<40);
T_2 = T(I_2);
P_2 = ch(I_2);

I_3 = find(T>46 & T<50);
T_3 = T(I_3);
P_3 = ch(I_3);

I_4 = find(T>50);
T_4 = T(I_4);
P_4 = ch(I_4);

f2 = figure;
f2.WindowState = 'maximized';

subplot(2,2,1)
plot(T_1,P_1)
grid on
xlabel('Time (s)')
ylabel('Amplitude')
title(['Signal on channel ',spec.channelnames{5},' - Partition ',num2str(1)])

subplot(2,2,2)
plot(T_2,P_2)
grid on
xlabel('Time (s)')
ylabel('Amplitude')
title(['Signal on channel ',spec.channelnames{5},' - Partition ',num2str(2)])

subplot(2,2,3)
plot(T_3,P_3)
grid on
xlabel('Time (s)')
ylabel('Amplitude')
title(['Signal on channel ',spec.channelnames{5},' - Partition ',num2str(3)])

subplot(2,2,4)
plot(T_4,P_4)
grid on
xlabel('Time (s)')
ylabel('Amplitude')
title(['Signal on channel ',spec.channelnames{5},' - Partition ',num2str(4)])


%% =============================================== Part 3: Another Channel
num = 17;
f3 = figure;
f3.WindowState = 'maximized';
ch_2  = sig(num,:);
plot(T,ch_2)
grid on
xlabel('Time (s)')
ylabel('Amplitude')
title(['Signal on channel ',spec.channelnames{num}])


%% ================================================== Part 4: Whole Signal
offset = max(max(abs(sig)))/2 ;
ElecName = spec.channelnames ;
disp_eeg(sig,offset,fs,ElecName) ;


%% ========================================================== Part 6: DFT
begs = [2,30,42,50];
len = 5;
times = cell(1,length(begs));
figs = cell(1,length(begs));
sigs = cell(1,length(begs));
dfts = cell(1,length(begs));

for i = 1:length(begs)
    idx = find(T>begs(i) & T<begs(i)+len);
    times{i} = T(idx);
    sigs{i} = ch(idx);
    
    per = 1/fs;
    N = length(idx); % time points
    dfts{i} = fft(sigs{i});
    
    Z2 = abs(dfts{i}/N);
    Z1 = Z2(1:floor(N/2)+1);
    Z1(2:end-1) = 2*Z1(2:end-1);
    
    figs{i} = figure;
    figs{i}.WindowState = 'maximized';
    
    subplot(2,1,1)
    plot(times{i},sigs{i})
    grid on
    xlabel('Time (s)')
    ylabel('Amplitude')
    title(['Signal on channel ',spec.channelnames{5},' - Partition ',num2str(i)])
    
    subplot(2,1,2)
    fr = fs*(0:floor(N/2))/N; % frequency range
    plot(fr,Z1) 
    grid on
    title(['Single-Sided Amplitude Spectrum of Partition ',num2str(i)])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
end


%% ======================================================== Part 7: Pwelch
f7 = figure;
f7.WindowState = 'maximized';
for i = 1:length(sigs)
    [pxx, fr] = pwelch(sigs{i},floor(fs/2),floor(fs/4),floor(fs/2),fs,'onesided');
    %length(fr)
    subplot(2,2,i)
    plot(fr,pxx) 
    grid on
    title(['Single-Sided Power Spectrum of Partition ',num2str(i)])
    xlabel('Frequency (Hz)')
    ylabel('Power')
end


%% ======================================================== Part 8: STFT
f8 = figure;
f8.WindowState = 'maximized';
window = 128;
noverlap = 64;
nfft = 128;
for i = 1:length(sigs)
    subplot(2,2,i)
    spectrogram(sigs{i},hamming(window),noverlap,nfft,fs,'onesided');
    title(['Single-Sided Spectogram of Partition ',num2str(i)])
end