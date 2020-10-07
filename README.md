MATLAB version 2014b (System: Intel Core i5-7200U @ 2.5GHz, 64-bit, 8 GB RAM) 

SIFT.m
- The input image can be changed in the image input section.
- The initial sections will show the Gaussian and DoG scale space. The DoG scale space is shown in negative for better visualisation.
- Next sections finds the keypoints and shows all keypoints mapped to the original image.
- Then, dominant orientation is found out and shown in the original image.
- Finally the descriptors are found out and are stored in the descriptor cell. 

Image Matching - Match_features.m
- The first section will read two images.
- The next section will find the descriptors using SIFT_func.m where the figures have been removed.
- Next section makes a matrix of all features in image 1 and another matrix of all features in image 2.
- Then features are matched and stored and the matched features are sorted based on the distances.
- In the next section 20 features with the minimum distances are plotted.
- There is another section which can be uncommented and the values of iteration variable(i) can be changed as given to show the 10 best matched features which were manually handpicked(based on visual correctness) for images 3_1.jpg,3_2.jpg and 5_1.jpg,5_2.jpg