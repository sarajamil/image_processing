% Comp Eng 4TN4
% Project 2
% Image Segmentation
% Sara Jamil
clear all; close all; clc;

% Input RBG image
% img = imread('backlit1.jpg');
img = imread('backlit2.jpg');
% img = imread('me-dark.jpg');

% Convert colour space
img_ntsc = rgb2ntsc(img);
orig_img = img_ntsc(:,:,1);

% % Original colour image
% figure;
% imshow(img);
% % Luminance  - gray-scale image
% figure;
% imshow(orig_img); 

[X,Y] = size(orig_img);
npixel = X*Y;

% Histogram
[hist, range] = imhist(orig_img);

% % Plot the histogram
% figure;
% imhist(orig_img)

%% Histogram Equalization

% Transfer function
p = hist./npixel;
Np = zeros(256,1);
for k = 1:256
    Np(k,1) = sum(p(1:k)); 
end
new_he = round(Np*255);

he_s = zeros(X,Y);

for i = 1:npixel
    [~, gray_level] = min(abs(range - orig_img(i)));
    he_s(i) = new_he(gray_level,1);
end

% figure;
% plot(imhist(HE_s/max(max(HE_s))));
% title('HE histogram');

% figure;
% imshow(uint8(HE_s));
% title('Histogram Equalization enhanced image - Gray-scale');

% Convert back to RGB output 
img_ntsc(:,:,1) = he_s / max(max(he_s));
out_he = ntsc2rgb(img_ntsc);

% figure;
% imshow(HE_out_img);
% title('Histogram Equalization enhanced image - Colour');


%% Transfer function S-curve for backlit images

% Estimate threshold using P1
T = 256/2;

x_data = [0,T-25:T+25,round((255-T)/2),round(255-(255-T)/3),255];
y_data = [0,T*ones(1,51),round((255-T)/2),round(255-(255-T)/3),255];
poly = polyfit(x_data,y_data,3);
new_trans = polyval(poly,0:255);

% % Transfer function
% figure;
% plot(trans_new_pixels);
% title('Transfer Function');
% xlabel('Input gray levels');
% ylabel('Output gray levels');

trans_s = zeros(X,Y);

for i = 1:npixel
    [~, gray_level] = min(abs(range - orig_img(i)));
    trans_s(i) = new_trans(1,gray_level);
end

% figure;
% plot(imhist(trans_s/max(max(trans_s))));
% title('S-curve histogram');

% figure;
% imshow(uint8(trans_s));
% title('Transfer function enhanced image - Gray-scale');

% Convert back to RGB output 
img_ntsc(:,:,1) = trans_s / max(max(trans_s));
out_trans = ntsc2rgb(img_ntsc);

% figure;
% imshow(trans_out_img);
% title('Transfer function enhanced image - Colour');

%% Backlit Image Restoration using OCTM

% Image segmentation
% Segment into under-exposed and over-exposed segments
% using optimal global thresholding
% Assume superposition of Gaussian distributions

GMModel = fitgmdist(hist,2); 
[mu_1,index1] = min(GMModel.mu);
[mu_2,index2] = max(GMModel.mu);
sigma_1 = GMModel.Sigma(:,:,index1);
sigma_2 = GMModel.Sigma(:,:,index2);
P1 = min(GMModel.PComponents);
P2 = max(GMModel.PComponents);

% Estimate T (using P1)
T = round(P1*256);

% Partion the original histogram
partition = ones(256,1);
partition(1:T,1) = 1;
partition(T+1:256,1) = 2;

p1 = zeros(256,1);
p2 = zeros(256,1);
p1_index = find(partition==1);
p2_index = find(partition==2);
p1(p1_index) = hist(p1_index);
p2(p2_index) = hist(p2_index);

over_exp = zeros(X,Y);
under_exp = zeros(X,Y);
L_mat = zeros(X,Y);

for i = 1:npixel
    [~, gray_level] = min(abs(range - orig_img(i)));
    L = P1*p1(gray_level) / (P1*p1(gray_level)+P2*p2(gray_level));
    L_mat(i) = L;
    if L <= 0.5
        over_exp(i) = gray_level-1;
    elseif L > 0.5
        under_exp(i) = gray_level-1;
    end
end

% Image of binary segmentation
% figure;
% imshow(L_mat); 
% title('Segmentation');
% 
% figure;
% imshow(uint8(over_exp));
% title('Over-exposed segment');
% 
% figure;
% imshow(uint8(under_exp));
% title('Under-exposed segment');


% Enhance over/under-exposed regions separately
% OCTM under-exposed image
output_under = OCTM(under_exp,p1,0.5,0.5);

% figure;
% imshow(uint8(output_under));
% title('OCTM of under-exposed segment');

% OCTM over-exposed image
output_over = OCTM(over_exp,p2,1.5,0.1);

% figure;
% imshow(uint8(output_over));
% title('OCTM of over-exposed segment');

% Linear combination of images
new_img = zeros(X,Y);
for i = 1:npixel
    new_img(i) = L_mat(i)*output_under(i)+((1-L_mat(i))*output_over(i));
end

% figure;
% imshow(uint8(new_img));
% title('OCTM enhanced image - Gray-scale');

% Convert back to RGB output
img_ntsc(:,:,1) = new_img / max(max(new_img));
out_octm = ntsc2rgb(img_ntsc);

% figure;
% imshow(out_octm);
% title('OCTM enhanced image - Colour');


%% Plots for comparison of backlit image restoration techniques

% Image Outputs
figure;  
subplot(2,2,1);
imshow(img);
title('Original');
subplot(2,2,2);
imshow(out_trans);
title('S-curve');
subplot(2,2,3);
imshow(out_he);
title('HE');
subplot(2,2,4);
imshow(out_octm);
title('OCTM');

% % Transfer functions
% X = reshape(orig_img*255,[npixel,1]);
% Y = reshape(uint8(new_img),[npixel,1]);
% 
% x = 0:255;
% y = 0:255;
% figure;
% plot(x,y,'k');
% hold on;
% plot(x,new_he,'r');
% hold on;
% plot(x,new_trans,'b');
% hold on;
% scatter(X,Y,'g');
% xlabel('Input Intensity');
% ylabel('Output Intensity');
% title('Transfer Function Comparison');
% legend('Original','HE','S-curve','OCTM');
% axis([0 255 0 255]);
% 

% % Histograms
% figure;
% subplot(2,2,1)
% plot(imhist(orig_img));
% title('Original');
% subplot(2,2,2)
% plot(imhist(trans_s/max(max(trans_s))));
% title('S-curve');
% subplot(2,2,3);
% plot(imhist(he_s/max(max(he_s))));
% title('HE');
% subplot(2,2,4);
% plot(imhist(new_img/max(max(new_img))));
% title('OCTM');


