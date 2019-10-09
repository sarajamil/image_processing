clear all; close all; clc;

% Import image
img = imread('C:\Users\Sara\Desktop\Comp Eng 4TN4\Pictures\owl-kitten.jpg');
figure(1)
imshow(img);
title('Original Image');

% Frequency domain of image
imgfft = fft2(double(img));
imgfft = fftshift(imgfft);
imgfftmag = abs(imgfft);
imgfftshow = mat2gray(log(imgfftmag+1));
figure(2)
imshow(imgfftshow);
title('2D FFT of Original Image');

% Meshgrid frequency matrices
% for computing freq domain filters
D0 = 32;    %cutoff frequency
D0 = D0^2;
[M,N,P] = size(img);
u = 0:(M-1);    %set up range of variables
v = 0:(N-1);

idx = find(u>(M/2));    %compute indices for meshgrid
idy = find(v>(N/2));
u(idx) = u(idx) - M;
v(idy) = v(idy) - N;

[U, V] = meshgrid(u,v); %meshgrid arrays

D = U.^2 + V.^2;    %compute distances to center of filter

% Ideal filter
idealH = double(D <= D0);
idealH = fftshift(idealH');
idealHshow = mat2gray(idealH);
figure(3)
imshow(idealHshow);
title('Ideal low-pass filter');

% Filtering the image with ideal LPF
idealF_img = zeros(size(imgfft));
for i = 1:P     %P=3 for colour images
    idealF_img(:,:,i) = idealH.*imgfft(:,:,i);
end
idealF_imgmag = abs(idealF_img);
idealF_imgshow = mat2gray(log(idealF_imgmag+1));
figure(4)
imshow(idealF_imgshow);
title('Frequency domain of ideal LP filtered image');

% Inverse FFT to obtain image in spatial domain
idealf_img = ifft2(ifftshift(idealF_img));
figure(5)
imshow(uint8(idealf_img));
title('Ideal LP filtered image');

% Gaussian filter
gaussH = exp(-(D.^2)./(2*(D0^2)));
gaussH = fftshift(gaussH');
gaussHshow = mat2gray(gaussH);
figure(6)
imshow(gaussHshow);
title('Gaussian low-pass filter');

% Filtering the image with gaussian LPF
gaussF_img = zeros(size(imgfft));
for j = 1:P
    gaussF_img(:,:,j) = gaussH.*imgfft(:,:,j);
end
gaussF_imgmag = abs(gaussF_img);
gaussF_imgshow = mat2gray(log(gaussF_imgmag+1));
figure(7)
imshow(gaussF_imgshow);
title('Frequency domain of gaussian LP filtered image');

% Inverse FFT to obtain image in spatial domain
gaussf_img = ifft2(ifftshift(gaussF_img));
figure(8)
imshow(uint8(gaussf_img));
title('Gaussian LP filtered image');
