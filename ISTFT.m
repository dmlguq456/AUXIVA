function [x] = ISTFT(X, Nfft, Nshift, len)

[nhfft, Nframe, Nch] = size(X);
assert(Nfft == (nhfft-1)*2,'Wrong FFT size')
    
Nol = Nfft-Nshift;
len_padd = Nframe*Nshift+Nol;
x = zeros(len_padd, Nch);

if Nshift == Nfft/4
    win = sqrt(2/3)*hanning(Nfft,'periodic');   %1/4shift
elseif Nshift == Nfft/2
    win = sin(pi*((0:1:Nfft-1)'+0.5)/Nfft);     %1/2 shift
end

for ch = 1 : Nch
    for frame = 1 : Nframe
        x_buff = real(ifft([X(:,frame,ch);conj(X(end-1:-1:2,frame,ch))]));
        x((frame-1)*Nshift+1:(frame-1)*Nshift+Nfft,ch) = x((frame-1)*Nshift+1:(frame-1)*Nshift+Nfft,ch) + x_buff.*win;       
    end
end
x = x(Nol+1:end,:);
x = x(1:len,:);
end