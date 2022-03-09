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
offset = max(max(abs(sig)))/2.5 ;
ElecName = spec.channelnames ;
disp_eeg(sig,offset,fs,ElecName) ;


%% =========================================================== Part 6: DFT
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
    xlim([0 fs/(2*3)])
end


%% ======================================================== Part 7: Pwelch
f7 = figure;
f7.WindowState = 'maximized';
window = floor(fs/4);
noverlap = floor(fs/8);
nfft = floor(fs);
for i = 1:length(sigs)
    subplot(2,2,i)
    [pxx, fr] = pwelch(sigs{i},window,noverlap,nfft,fs,'onesided');
    plot(fr,pxx) 
%     pspectrum(sigs{i},fs)
    grid on
    title(['Single-Sided Power Spectrum of Partition ',num2str(i)])
    xlabel('Frequency (Hz)')
    ylabel('Power')
end


%% ========================================================== Part 8: STFT
f8 = figure;
f8.WindowState = 'maximized';
window = 128;
noverlap = 64;
nfft = 128;
for i = 1:length(sigs)
    subplot(2,2,i)
    spectrogram(sigs{i},hamming(window),noverlap,nfft,fs,'onesided');
    xlim([0 fs/(2*Down)])
    %yticks(linspace(begs(i),begs(i)+5,10))
    title(['Single-Sided Spectogram of Partition ',num2str(i)])
end


%% ================================================== Part 9: Downsampling
part    = 2; % partition number
Down    = 3;
z_lpf   = lowpass(sigs{part},fs/(2*Down),fs);
z_down  = downsample(z_lpf,Down);

time = downsample(times{part},Down);
N = length(z_down); % time points
Y = fft(z_down);

Z2 = abs(Y/N);
Z1 = Z2(1:floor(N/2)+1);
Z1(2:end-1) = 2*Z1(2:end-1);

f91 = figure;
f91.WindowState = 'maximized';

subplot(2,1,1)
plot(time,z_down)
grid on
xlabel('Time (s)')
ylabel('Amplitude')
title(['Downsampled Signal on channel ',spec.channelnames{5},' - Partition ',num2str(part)])

subplot(2,1,2)
fr = (fs/Down)*(0:floor(N/2))/N; % frequency range
plot(fr,Z1) 
grid on
title(['Single-Sided Amplitude Spectrum of Downsampled Partition ',num2str(part)])
xlabel('Frequency (Hz)')
ylabel('Amplitude')

f92 = figure;
f92.WindowState = 'maximized';
window = floor(128/Down);
noverlap = floor(64/Down);
nfft = floor(128/Down);
spectrogram(z_down,hamming(window),noverlap,nfft,fs/Down,'onesided');
title(['Single-Sided Spectogram of Downsampled Partition ',num2str(part)])