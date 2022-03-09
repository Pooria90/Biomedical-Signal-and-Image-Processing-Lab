clear; clc;

%% ===================================================================== 1
sig = load('X_org.mat');
sig = sig.X_org;
load('Electrodes') ;
offset = max(abs(sig(:)));
fs = 250;
ElecName = Electrodes.labels ;
disp_eeg(sig,offset,fs,ElecName);
title('Original Signal')


%% ===================================================================== 2
sig_n=load('X_noise.mat');
sig_n=sig_n.X_noise;
load('Electrodes') ;
offset = max(abs(sig_n(:)));
disp_eeg(sig_n,offset,fs,ElecName);
title('Noisy Signal')


%% ===================================================================== 3
