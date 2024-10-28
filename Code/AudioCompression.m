%%
[A,fs]=audioread('sample.wav');
%% Sender
C=[];
A=A';
for i=512:512:numel(A)
    B=DCT(A(i-511:i));
    C=[C, B(1:128)];
end
%% Reciever
A2=[];
for i=128:128:numel(C)
    S=[C(i-127:i),zeros(1,384)];
    S=IDCT(S); 
    A2=[A2,S];
end
%% Evaluation
dis=numel(A)-numel(A2);
A2=[A2,zeros(1,dis)];
PSNR=psnr(A2,A);
MSError=immse(A2,A);
SNR=snr(A2,A);
fprintf("PSNR %f\n", PSNR)
fprintf("MSE Error %f\n", MSError)
fprintf("SNR %f\n", SNR)
%% plot
figure, 
title('Audio Compression By DCT');
subplot(3,1,1)
plot(A);
title('Original Audio Signal')
xlabel('Samples');
ylabel('Amp.');

subplot(3,1,2)
plot(A2,'r')
title('Compressed Audio Signal')
xlabel('Samples');
ylabel('Amp.');

subplot(3,1,3)
plot((A2-A),'y')
title("Difference in original and compressed signal")
xlabel('Samples');
ylabel('Amp.');

%% dct1d function
function b=DCT(a)
temp=a;
if (size(temp,1) == 1)
    a = a(:); 
end
[n,m] =size(a);
aa = a(1:n,:);
ww = (exp(-1i*(0:n-1)*pi/(2*n))/sqrt(2*n)).';
ww(1) = ww(1) / sqrt(2);
if rem(n,2)==1 || ~isreal(a)
    % intermediate matrix
    y = zeros(2*n,m);
    y(1:n,:) = aa;
    y(n+1:2*n,:) = flipud(aa);
    yy = fft(y);
    yy = yy(1:n,:);
else
    y = [ aa(1:2:n,:); aa(n:-2:2,:) ];
    yy = fft(y);
    ww = 2*ww;  % Double the weights
end
b = ww(:,ones(1,m)) .* yy;
if isreal(a)
    b = real(b); 
end
if (size(temp,1) == 1) 
    b = b.'; 
end

end

%% IDCT
function a = IDCT(b)
 if size(b,2)>1
        do_trans = 1;
    else
        do_trans = 0;
    end
b = b(:);
[n,m]= size(b);
 bb = b(1:n,:);
if rem(n,2)==1 || ~isreal(b)
    ww = sqrt(2*n) * exp(1i*(0:n-1)*pi/(2*n)).';
    ww(1) = ww(1) * sqrt(2);
    W = ww(:,ones(1,m));
    yy = zeros(2*n,m);
    yy(1:n,:) = W.*bb;
    yy(n+2:n+n,:) = -1i*W(2:n,:).*flipud(bb(2:n,:));
    y = ifft(yy);
    a = y(1:n,:);

else 
    % Compute precorrection factor
    ww = sqrt(2*n) * exp(1i*pi*(0:n-1)/(2*n)).';
    ww(1) = ww(1)/sqrt(2);
    W = ww(:,ones(1,m));
    y = ifft(W.*bb);
    a = zeros(n,m);
    a(1:2:n,:) = y(1:n/2,:);
    a(2:2:n,:) = y(n:-1:n/2+1,:);
end
if isreal(b)
    a = real(a); 
end
if do_trans
    a = a.'; 
end
end