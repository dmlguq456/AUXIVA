clc; clear;

%% data read
[x, fs] = audioread('./input.wav');

%% parameter
nfft = 2048;
iteration = 50;

%% Processing
[y] = AUXIVA(x,nfft,iteration);

%% Spectrogram
figure('Position',[100 100 1000 300]);
spectrogram(x(:,1),512,384,512,fs,'yaxis'); colormap jet; caxis([-140  -30]); xlabel('Time(s)'); ylabel('Frequency(kHz)');

figure('Position',[100 100 1000 700]);
subplot(2,1,1);
spectrogram(y(:,1),512,384,512,fs,'yaxis'); colormap jet; caxis([-140  -30]); xlabel('Time(s)'); ylabel('Frequency(kHz)');
subplot(2,1,2);
spectrogram(y(:,2),512,384,512,fs,'yaxis'); colormap jet; caxis([-140  -30]); xlabel('Time(s)'); ylabel('Frequency(kHz)');

audiowrite('output1.wav',y(:,1),fs);
audiowrite('output2.wav',y(:,2),fs);