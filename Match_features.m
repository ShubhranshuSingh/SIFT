%% (e) Image Matching 
% Takes at least 5min to run
clc
clear
%% Read Images
I1 = imread('images/3_1.jpg');
I2 = imread('images/3_2.jpg');

%% Find the descriptors of the 2 image
D1 = SIFT_func(I1);
D2 = SIFT_func(I2);

%% Save features from all octaves into one matrix
features1=[];
for i=1:size(D1,1)
    if(size(D1{i})~=0)
        % Coordinates multiplied by 2^(octave-1) to find their location in
        % the Original image
        features1 = [features1; [2^(i-1)*D1{i}(:,1:2) D1{i}(:,3:end)]];
    end
end
features2=[];
for i=1:size(D2,1)
    if(size(D2{i})~=0)
        % Coordinates multiplied by 2^(octave-1) to find their location in
        % the Original Image
        features2 = [features2; [2^(i-1)*D2{i}(:,1:2) D2{i}(:,3:end)]];
    end
end

%% Find all matched points (Takes time to run as every feature in 1st image is checked with every feature in 2nd image)
% This matrix contains y1,x1,y2,x2 and dist of the matched descriptors 
matched_points = [];
for i=1:size(features1,1)
    min_ = Inf;
    for j=1:size(features2,1)
        dist = Norm2(features1(i,3:end)-features2(j,3:end)); % Distance measure is norm 2.
        if(dist<min_)
            matched_point = [features1(i,1:2),features2(j,1:2)];
            min_ = dist;
        end
    end
    matched_points = [matched_points; [matched_point, min_]];
end
% Sort with respect to the distance of the two features
[~,idx] = sort(matched_points(:,5));
matched_points = matched_points(idx,:);

%% Display all 20 best matched points which have minimum distance


showMatchedFeatures(imresize(I1,[NaN 512]),imresize(I2,[NaN 512]),[matched_points(1:20,2),matched_points(1:20,1)],[matched_points(1:20,4),matched_points(1:20,3)],'montage','Parent',axes);

%% Display 10 best matched which are handpicked
% If I1 is 3_1.jpg and I2 is 3_2.jpg use i= [1,2,3,5,8,9,10,11,12,13] instead of
% i=1:1
% If I1 is 5_1.jpg and I2 is 5_2.jpg use i= [2,4,12,18,22,24,25,28,42,47] instead of
% i=1:1
% Uncomment the code below
% for i=1:1
%     figure(); 
%     showMatchedFeatures(imresize(I1,[NaN 512]),imresize(I2,[NaN 512]),[matched_points(i,2),matched_points(i,1)],[matched_points(i,4),matched_points(i,3)],'montage','Parent',axes);
% end