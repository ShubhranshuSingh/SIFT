clc
clear

% Constants
n = 4; % Number of octaves
s = 5; % Number of scales. Should be greater than 3 
k = sqrt(2); % Scaling factor
sigma = 1/sqrt(2); % Standard Deviation


%% (a) Scale Space Decomposition

%% Image input
I = cell(n,1);
I{1} = imread('images/3_1.jpg');
I{1} = double(rgb2gray(imresize(I{1},[NaN 512])));

% Resizing image for different octaves
for i=2:n
    I{i} = imresize(I{i-1},0.5);
end

% Value of standard deviation at different octaves
sigma_octave = zeros(n,s);

for i=1:n
    for j=1:s
        sigma_octave(i,j) = sigma*k^(2*(i-1)+j-1);
    end
end

% Sizes of different octaves
size_octave = zeros(n,2);
for i=1:n
    size_octave(i,:) = size(I{i});
end

%% Gaussian Scale Space
L = cell(n, s);

for i=1:n
    for j=1:s
        % Size of filter varied with standard deviation
        L{i,j} = convolution_handmade(I{i},Gaussian(sigma_octave(i,j),ceil(6*sigma_octave(i,j))));
    end
end

%% DoG Scale Space

D = cell(n,s-1);

% D_pad to be used to detect keypoints by mirror padding
D_pad = cell(n,s-1);

for i=1:n
    for j=1:s-1
        D{i,j} = L{i,j+1}-L{i,j};
        D_pad{i,j} = padreflect(D{i,j});
    end
end

%% Show Scale Space Images
figure;
count = 1;
for i=1:n
    for j=1:s
        subplot(n,s,count),imshow(uint8(L{i,j}))
        title(strcat(strcat('Octave:',int2str(i)),strcat(' Scale:',int2str(j))));
        count = count+1;
    end
end
suptitle('Gaussian Blurred Images');

figure;
count = 1;
for i=1:n
    for j=1:s-1
        % Displaying negative of the DoG image for better visualisation
        subplot(n,s-1,count),imshow(255-uint8(D{i,j}))
        title(strcat(strcat('Octave:',int2str(i)),strcat(' Scale:',int2str(j))));
        count = count+1;
    end
end
suptitle('DoG Images');
%% (b) Keypoint Detection

K = cell(n,s-3);  

for i=1:n
    for j=1:size_octave(i,1)
        for k=1:size_octave(i,2)
            for l=2:s-2
                % Finding a 3x3x3 volume around each point
                temp =  [D_pad{i,l-1}(j:j+2,k:k+2) D_pad{i,l}(j:j+2,k:k+2) D_pad{i,l+1}(j:j+2,k:k+2)];
                [~,~,min_temp] = minn(temp); %minimum value in the volume
                [~,~,max_temp] = maxx(temp); %maximum value in the volume

                % Storing if found to be the keypoint
                if(D{i,l}(j,k) == min_temp || D{i,l}(j,k) == max_temp)
                    K{i,l-1} = [K{i,l-1}; [j k]];
                end
            end
        end
    end
end

%% Show all keypoints on the original image

% Different colors for different octave of the keypoint
% The following function was found online. Link- https://in.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
colors = distinguishable_colors(n,{'w','k'});
figure;
imshow(uint8(I{1}));
hold on;
for i=1:n
    y = [];
    x = [];
    for j=1:s-3
        if(size(K{i,j})~=0)
            % Coordinates multiplied by 2^(octave-1) to find their location
            % in the Original Image
            y = [y; (K{i,j}(:,1))*2^(i-1)]; 
            x = [x; (K{i,j}(:,2))*2^(i-1)];
        end
    end
    plot(x,y,'.','Color', colors(i,:),'DisplayName',strcat('Octave:',int2str(i)));
end
legend('show')

%% (c)Orientation assignment- Step 1
% M and Theta calculated for all points in Gaussian Scale Space

M = cell(n,s-3);
Theta = cell(n,s-3);

for i=1:n
    for j=1:s-3
        M{i,j} = zeros(size_octave(i,1),size_octave(i,2));
        Theta{i,j} = zeros(size_octave(i,1),size_octave(i,2));
    end
end

for i =1:n
    for j=2:size_octave(i,1)-1
        for k=2:size_octave(i,2)-1
            for l=3:s-1
                % y derivative
                y_ = L{i,l}(j+1,k)-L{i,l}(j-1,k);
                % x derivative
                x_ = L{i,l}(j,k+1)-L{i,l}(j,k-1);
                % Store the magnitude of derivative
                M{i,l-2}(j,k) = ((x_)^2+(y_)^2)^0.5;
                if(x_ == 0 && y_ == 0)
                    theta = 0;
                elseif(x_>0 && y_>0)
                    theta = atand(y_/x_);
                elseif(x_>0 && y_<0)
                    theta = 360-atand(abs(y_/x_));
                else
                    theta = 180 + atand(y_/x_);
                end
                % Store the theta
                Theta{i,l-2}(j,k) = theta;
            end
        end
    end
end

%% Orientation assignment- Step 2
O = cell(n,s-3);

for i =1:n
    for l=2:s-2
        xy = K{i,l-1};
        for j=1:size(xy)
            % To check if a 16x16 neighborhood exists
            % If not then the keypoint is removed
            if(xy(j,1)>=8 && xy(j,1)<=size_octave(i,1)-8 && xy(j,2)>=8 && xy(j,2)<=size_octave(i,2)-8)
                % Histogram
                bins = zeros(36,1);
                gaussian_weight = Gaussian16(1.5*sigma_octave(i,l+1));
                % Looking in a 16x16 neighborhood
                for a=-7:8
                    for b=-7:8
                        theta = Theta{i,l-1}(a+xy(j,1),b+xy(j,2));
                        % Storing M times gaussian weight in a correct bin
                        if(ceil(theta/10)==0)
                            bins(1) = bins(1)+M{i,l-1}(a+xy(j,1),b+xy(j,2))*gaussian_weight(a+8,b+8);
                        else
                            bins(ceil(theta/10)) = bins(ceil(theta/10))+ M{i,l-1}(a+xy(j,1),b+xy(j,2))*gaussian_weight(a+8,b+8);
                        end
                    end
                end
                % Storing location of keypoint with its dominant orientation
                [theta_dom,~,~] = maxx(bins);
                theta_dom = theta_dom*10;
                O{i,l-1} = [O{i,l-1}; [xy(j,1) xy(j,2) theta_dom]];
            end
        end
    end
end

%% Display keypoints with orientations

% Magnitude of vector indicates the octave it was found at
% Longer vector implies larger octave
figure;
imshow(uint8(I{1}));

for i=1:n
    x = [];
    y = [];
    delx = [];
    dely = [];
    for l=2:s-2
          if(size(O{i,l-1},1)~=0)
              % Coordinates multiplied by 2^(octave-1) to find their
              % location in the Original Image
              x = [x; (O{i,l-1}(:,2))*2^(i-1)];
              y = [y; (O{i,l-1}(:,1))*2^(i-1)];
              delx = [delx; ((i+1)^3)*cosd(O{i,l-1}(:,3))];
              dely = [dely; ((i+1)^3)*sind(O{i,l-1}(:,3))];
          end
    end
    hold on;
    plot(x,y,'.','Color', colors(i,:),'DisplayName',strcat('Octave:',int2str(i)));
    hold on;
    quiver(x,y,delx,dely,0,'Color', colors(i,:),'DisplayName',strcat('Octave:',int2str(i)));
end
legend('show');



%% (d) SIFT 128D descriptor

descriptor = cell(n,1);

for i =1:n
    for l=2:s-2
        xy = O{i,l-1};
        for j=1:size(xy)
            % Bins for each 4x4 patch
            bins = cell(4,4);
            for a=1:4
                for b=1:4
                    bins{a,b} = zeros(1,8);
                end
            end
            gaussian_weight = Gaussian16(8);

            % 16x16 neigborhood
            for a=-7:8
                for b=-7:8
                    theta = Theta{i,l-1}(a+xy(j,1),b+xy(j,2));
                    % Storing M times gaussian weight in a correct bin
                    if(ceil(theta/45)==0)
                        bins{ceil(a/4)+2,ceil(b/4)+2}(1) = bins{ceil(a/4)+2,ceil(b/4)+2}(1)+M{i,l-1}(a+xy(j,1),b+xy(j,2))*gaussian_weight(a+8,b+8);
                    else
                        bins{ceil(a/4)+2,ceil(b/4)+2}(ceil(theta/45)) = bins{ceil(a/4)+2,ceil(b/4)+2}(ceil(theta/45))+ M{i,l-1}(a+xy(j,1),b+xy(j,2))*gaussian_weight(a+8,b+8);
                    end
                end
            end
            vect = [];
            
            % Making a 128D vector
            for a=1:4
                for b=1:4
                    vect = [vect, bins{a,b}];
                end
            end

            % Normalise the 128D vector
            norm_vect = 0;
            for a=1:128
                norm_vect = norm_vect+vect(a)^2;
            end
            vect = vect./sqrt(norm_vect);

            % Threshold so that max value is 0.2
            for a=1:128
                if(vect(a)>0.2)
                    vect(a)=0.2;
                end
            end

            % Renormalise
            norm_vect = 0;
            for a=1:128
                norm_vect = norm_vect+vect(a)^2;
            end
            vect = vect./sqrt(norm_vect);

            % Storing x,y and 128D vector at appropriate octave
            descriptor{i} = [descriptor{i};[xy(j,1) xy(j,2)], vect];

        end
    end
end