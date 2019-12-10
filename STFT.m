function [X] = STFT(x, Nfft, Nshift)

[wavlen, Nch] = size(x);
Nframe = ceil((wavlen)/Nshift);
len_padd = Nframe*Nshift;
x = [x;zeros(len_padd-wavlen,Nch)];

X = zeros(Nfft/2+1,Nframe,Nch);

if Nshift == Nfft/4
    win = sqrt(2/3)*coder.const(@hanning,Nfft,'periodic');%1/4shift
elseif Nshift == Nfft/2
    win = sin(pi*((0:1:Nfft-1)'+0.5)/Nfft); %1/2 shift
end

for ch = 1 : Nch
    x_buff = zeros(Nfft,1);
    for frame = 1 : Nframe
        x_buff = [x_buff(Nshift+1:end);x((frame-1)*Nshift+1:frame*Nshift,ch)];
        X_temp = fft(x_buff.*win);
        X(:,frame,ch) = X_temp(1:Nfft/2+1);
    end
end

end