function kernel = Gaussian16(sigma)
    %% Returns a Gaussian kernel with size 16x16
    % Intialise with zeros
    kernel = zeros([16 16]);
    % Sampling 2D Gaussian Function at discrete points 
    for i = -7:8
        for j = -7:8
            kernel(i+8,j+8) = exp(-(i^2+j^2)/(2*(sigma^2)));
        end
    end
    
    % Normalisation
    kernel = kernel./summ(kernel);
end