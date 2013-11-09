% Article - Encodage Spatial basé sur les propriétés du
%           système visuel humain. Jérôme Schmaltz, MTI810.
%
% Utilisation des pyramides gaussiennes pour simuler un
% comportement fovéal selon un point de fixation.
%
% Parameters:
%
% x         : foveation point x coordinate
% y         : foveation point y coordinate
% imagepath : image source path
%
% return the processed image.
%
% example : gpyramid('juliette_pommes.jpg', 50, 50);
%
function [ image ] = gpyramid( imagepath, x, y )

% read the image
image = double(imread(imagepath));

% determine the image size
[size_x, size_y] = size(image(:,:,1));

% interpolate the image between 0-1.
max_cval = 1.0;
range = [min(image(:)) max(image(:))];
image = max_cval.*(image-range(1)) ./ (range(2)-range(1));

image_ycbcr = rgb2ycbcr(image);

% spatial frequency decay constant
% based on HVS best fir parameter
alpha = 0.106;

% half-resolution eccentricity constant
e2 = 2.3;

% minimal contrast threshold
CT0 = 1 / 64;

% viewing distance
% between retina and image plane
% in meters
V = 0.3;

% monitor dot pitch
% in meter
D = 0.25 * 10^(-3);

% distance point-foveration point
% (will be calculate for each pixels)
% in meters
l = 0;

% eccentricity at point c(x,y)
% calculated for each pixel of the image
ec = 0;

% pyramid level required to be sent at each 
% point of the image
pyrlevel = 0;
pyrmatrix = zeros(size_x, size_y);

% loop through the image pixels
for i = 1 : size_x
    for j = 1 : size_y

        % calculate distance foveation point - pixel
        l = ( (i-x)^2 + (j-y)^2 )^(1/2);

        % calculate the eccentricity
        ec = 180 * atan(D * l / V) / pi;

        % display resolution
        fm = pi/((atan(D*(l+1)/V)-atan(D*(l-1)/V))*180);        

        % calculate the cutoff frequency
        fc = (e2 * log(1/CT0)) / (alpha * (ec + e2));

        % calculate the pyramide level
        pyrlevel = fm / fc;

        % store the pyramidal level for
        % each pixel
        pyramid_level_mat(i, j) = pyrlevel;

    end
end

% translate the pyramidal level between [1 to (max_level-1)]
min_level = min(min(pyramid_level_mat));
max_level =  (max(max(pyramid_level_mat)) - min_level) + 1;

pyramid_level_mat_res = 1 - (pyramid_level_mat ./ (max_level - min_level));
pyramid_level_mat = (pyramid_level_mat - min_level) + 1;

% display the pyramid level graph
figure;
surfc(pyramid_level_mat);
shading interp
colormap('jet');

% display the distribution resolution 
figure;
surfc(pyramid_level_mat_res);
shading interp
colormap('gray');

max_pyr_level = ceil(max_level);
pyramid_images = ones(size_x, size_y, 3, max_pyr_level);

pyramid_images(:, :, 1, 1) = image_ycbcr(:, :, 1);
pyramid_images(:, :, 2, 1) = image_ycbcr(:, :, 2);
pyramid_images(:, :, 3, 1) = image_ycbcr(:, :, 3);

for comp=1:3
    reduced_image = image_ycbcr(:, :, comp);
    for k = 2:max_pyr_level
        % image is reduced
        reduced_image = impyramid(reduced_image, 'reduce');
        % and is expanded to be the same
        % dimensions as the original image (linear interpolation)
        level_image = imresize(reduced_image, 2^(k-1), 'nearest', 'OutputSize', [size_x, size_y]);
        % store the image
        pyramid_images(:, :, comp, k) = level_image;
    end
end

% creating the final image from multiresolution
% pyramidal leveled images.
final_image_ycbcr = zeros(size_x, size_y, 3);

for comp=1:3
    for l=1:size_x
        for m=1:size_y
            
            % get the discrited level
            level_up = ceil(pyramid_level_mat(l, m));
            level_low = floor(pyramid_level_mat(l, m));
            
            % interpolate the value

            % get the related image
            level_image_up = pyramid_images(:, :, comp, level_up);
            level_image_low = pyramid_images(:, :, comp, level_low);
            
            color_up = level_image_up(l, m);
            color_low = level_image_low(l, m);
            
            final_color =  color_up + ((color_up - color_low) / (level_up - level_low) * (pyramid_level_mat(l, m) - level_low));
            
            % get the related pixel
            final_image_ycbcr(l, m, comp) = final_color;
        end 
    end
end

% show the final image
final_image_rgb = ycbcr2rgb(final_image_ycbcr);

% readjust the range of the final image in order to ensure values
% 0.0 < foveated(i,j) < 1.0
range = [min(final_image_rgb(:)) max(final_image_rgb(:))];
final_image_rgb = (final_image_rgb-range(1)) ./ (range(2) - range(1));

figure;
imshow(final_image_rgb);

% return the final image.
gpyramid = final_image_rgb;

end

