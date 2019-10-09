clc; clear all; close all; 

inputfile = 'fancycar';
inputext = '.jpg';
input = [inputfile,inputext];
input_img = imread(input);
[rows, columns, numberOfColorChannels] = size(input_img);

if numberOfColorChannels == 3
    input_img = rgb2gray(input_img);
end

array_hist = imhist(input_img);   % probability from the histogram

% constraint parameters
L = 256;    
Lbar = 256;
[x,y] = size(input_img);    %size of the input image is extracted
npixels = x*y;         %the total number of pixels calculated
probability = array_hist/npixels;  %probability of the occurence of a pixel
diff = 2; % for average intensity constraint

% constraints
I=[1:L];

A = [ones(1,256)                        % fundamental constraints
    triu(repmat(probability',L,1))];	% average intensity constraint

b = [255                                % fundamental constraints
    diff+probability.*I'];              % average intensity constraint

% linprog arguements
Aeq = []; % empty matrix
Beq = []; % empty matrix
distortion = 3;    % distortion factor
lower_bound = (1/distortion)*ones(1,L); 
upper_bound = Lbar*ones(1,L);
% calculation of the step increment
Sj = linprog(-probability,A,b,Aeq,Beq,lower_bound,upper_bound);
output = input_img; % input image is set up as a new variable in order to be used in the re allocation of pixels.
% for loop, used to re-allocate the pixels
for i=1:x
    for j=1:y
        output(i,j) = (sum(Sj(1:input_img(i,j)))+0.5);
    end
end

name = ['octm_avgint_',inputfile,'_d',num2str(distortion),'_diff',num2str(diff)];
figure(1)  
subplot(2,2,1);
imshow(input_img);
title('Original Image');
subplot(2,2,2)
plot(imhist(input_img));
title('Original Histogram');
subplot(2,2,3);
imshow(output);
title('OCTM Output');
subplot(2,2,4);
plot(imhist(output));
title('Output Histogram');
% print(name,'-dpng');

% TRANSFER FUNCTION
% real one using Sj values
X=(1:L)';
Y=zeros(L,1);
for i=1:256
    Y(i,1)=sum(Sj(1:i));
end

figure(2)
plot(X,Y);
title('Transfer function');
xlabel('Input gray levels');
ylabel('Output gray levels');
% print([name,'_TRANS2'],'-dpng');
