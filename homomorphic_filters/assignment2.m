% 4TN4 Assignment 2
% Homomorphic filter
clear all; close all; clc;

% img = imread('lena.png');
% img = imread('street_orig.png');
% img = imread('girl_orig.png');
img = imread('cave.jpg');
imsize = size(img);
if (numel(imsize)>2)
    img = rgb2gray(img);
end

figure(1)
imshow(img);
title('Original Image');

cimg = im2double(img);
%add 1 to pixels to remove 0 values resulting in undefined log values
cimg = cimg + 1;

% Natural logarithm
limg = log(cimg);
% figure(3)
% imshow(limg);
% title('Natural Logarithm');

% Fourier Transform
fimg = fft2(limg);
% figure(4)
% imshow(uint8(abs(fftshift(fimg))));
% title('Fourier Transform');

% Gaussian High-Pass Filter
D0 = 128;    %cutoff frequency
D0 = D0^2;
[M,N,P] = size(img);
u = 0:(M-1);    %set up range of variables
v = 0:(N-1);

idx = find(u>(M/2));    %compute indices for meshgrid
idy = find(v>(N/2));
u(idx) = u(idx) - M;
v(idy) = v(idy) - N;

%meshgrid frequency matrices for computing freq domain filters
[U, V] = meshgrid(u,v); 

D = U.^2 + V.^2;    %compute distances to center of filter

gaussH = 1 - (exp(-(D.^2)./(2*(D0^2)))); %gaussian lowpass filter

% Homomorphic Filter
gammaL = 0.9;
gammaH = 5.5;
diff = gammaH - gammaL;

gaussH = diff*gaussH + gammaL;
% figure(5)
% imshow(mat2gray(fftshift(gaussH)));
% title('Filter');

himg = gaussH'.*fimg;
% figure(6)
% imshow(uint8(abs(fftshift(himg))));
% title('Filtered Image - Freq Domain');

% Inverse Fourier Transform
ifimg = ifft2(himg);
% figure(7)
% imshow(ifimg);
% title('Inverse Fourier Tranform');

% Exponent
eimg = exp(ifimg) - 1;

figure(2)
imshow(eimg);
title('Homomorphic Filtering');