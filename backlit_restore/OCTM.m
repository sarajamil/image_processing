function output_img = OCTM(orig_img,hist,distortion,diff)
% Optimal Contrast Tone Mapping

[X,Y,Z] = size(orig_img);
npixel = X*Y;

% Parameters   
L = 256;
Lbar = 256;
probability = hist/npixel;  %proabbility of the occurance of a pixel

% Constraints
I = [1:L];

A = [ones(1,L)                                  % fundamental constraints
    triu(repmat(probability',L,1))];            % average intensity constraint

b = [L-1                                        % fundamental constraints
    diff+probability.*I'];                      % average intensity constraint

Aeq = [];
Beq = [];
lower_bound = (1/distortion)*ones(1,L); 
upper_bound = Lbar*ones(1,L);

Sj = linprog(-probability,A,b,Aeq,Beq,lower_bound,upper_bound);
output_img = orig_img; 

for i=1:X
    for j=1:Y
        output_img(i,j) = sum(Sj(1:orig_img(i,j))+0.5);
    end
end

end