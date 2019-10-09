clear all; close all; clc;

inputfile = 'low-contrast';
inputext = '.png';
input = [inputfile,inputext];
InputImage = imread(input);

[rows, columns, numberOfColorChannels] = size(InputImage);
 
% Convert to grayscale image
if numberOfColorChannels == 3
    InputImage = rgb2gray(InputImage); 
end

pixelTotal=size(InputImage,1)*size(InputImage,2);

% Calculating number of pixels for each occurence of a gray level
freq=zeros(256,1);

for i=1:size(InputImage,1)
    for j=1:size(InputImage,2)
        value=InputImage(i,j);
        freq(value+1)=freq(value+1)+1;
    end
end

sum=0;
no_bins=255;

% Cumulative distribution probability
cdf=zeros(256,1);
s=zeros(256,1);

for i=1:256
   sum=sum+freq(i);
   cdf(i)=sum;
   s(i)=round((cdf(i)/pixelTotal)*no_bins); %transfer function
end

equalizedImage=uint8(zeros(size(InputImage,1),size(InputImage,2)));

for i=1:size(InputImage,1)
    for j=1:size(InputImage,2)
            equalizedImage(i,j)=s(InputImage(i,j)+1); 
            %mapping gray level of original image to equalized image
    end
end

name = ['histeq_',inputfile];

figure('Position', [100, 425, 1400, 400])
subplot(1,2,1);
imshow(InputImage);
title('Original Image');
subplot(1,2,2);
plot(imhist(InputImage));
title('Histogram of Original');
print([name,'_original'],'-dpng');
figure('Position', [100, 0, 1400, 400])
subplot(1,2,1);
imshow(equalizedImage);
title('Histogram equalization');
subplot(1,2,2);
plot(imhist(equalizedImage));
title('Histogram of Equalized Image');
print([name,'_enhanced'],'-dpng');
