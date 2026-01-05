clear; clc; close all;
fileName = 'TennisSet1.zip';
zipfileName = 'SET1';
unzip(fileName, zipfileName);
photo = imread(fullfile(zipfileName, 'stennis.0.ppm'));

grayPhoto = rgb2gray(photo);
grayDouble = im2double(grayPhoto);
photoDouble = im2double(photo); 

grayFiltered = medfilt2(grayPhoto, [3 3]);

whiteDet = (photoDouble(:,:,1) > 0.8) & (photoDouble(:,:,2) > 0.8) & (photoDouble(:,:,3) > 0.8);
lowDet = imgaussfilt(grayFiltered, 5);
detailsCont = (grayFiltered - lowDet) > 0.1;

ballMaskR = whiteDet & detailsCont;
ballMaskR = imopen(ballMaskR, strel('disk',2));
ballMaskR = imclose(ballMaskR, strel('disk',2));
ballMaskR = bwareaopen(ballMaskR, 10);  

[balls, numBalls] = bwlabel(ballMaskR);

stats = regionprops(balls, 'Area', 'Centroid', 'BoundingBox');
if ~isempty(stats)
    [~, idx] = max([stats.Area]); 
    ballMask = (balls == idx);
    fprintf('Ball detected!\n');
else
    ballMask = false(size(grayPhoto));
    fprintf('No ball detected.\n');
end

redMask = (photoDouble(:,:,1) > 0.6) & (photoDouble(:,:,2) < 0.4) & (photoDouble(:,:,3) < 0.4);
blackMask = (photoDouble(:,:,1) < 0.35) & (photoDouble(:,:,2) < 0.35) & (photoDouble(:,:,3) < 0.35);

batMaskR = redMask | blackMask;

batMaskR = imopen(batMaskR, strel('rectangle',[3 3]));
batMaskR = imclose(batMaskR, strel('rectangle',[3 3]));
batMaskR = bwareaopen(batMaskR, 50);

[bats, numBats] = bwlabel(batMaskR);
batStats = regionprops(bats, 'Area', 'MajorAxisLength', 'MinorAxisLength');

batCand = false(size(batMaskR));
minBatArea = 200; maxBatArea = 2000;
minBat = 3;

for k = 1:numBats
    ratio = batStats(k).MajorAxisLength / batStats(k).MinorAxisLength;
    if batStats(k).Area >= minBatArea && batStats(k).Area <= maxBatArea && ratio >= minBat
        batCand(bats == k) = true;
    end
end

[final, ~] = bwlabel(batCand);
propsF = regionprops(final, 'Area');

if ~isempty(propsF)
    [~, idx] = max([propsF.Area]);   
    batMask = (final == idx);
else
    batMask = false(size(batMaskR));
end

if any(ballMask(:))

    ballProps = regionprops(ballMask, ...
        'Area', ...
        'Centroid', ...
        'BoundingBox', ...
        'MajorAxisLength', ...
        'MinorAxisLength', ...
        'Perimeter', ...
        'Extrema', ...
        'Eccentricity');   

    ballArea     = ballProps.Area;
    ballCentroid = ballProps.Centroid;      
    ballBBox     = ballProps.BoundingBox;     
    ballMajor    = ballProps.MajorAxisLength;
    ballMinor    = ballProps.MinorAxisLength;
    ballPerim    = ballProps.Perimeter;
    ballExtrema  = ballProps.Extrema;         
    ballEcc      = ballProps.Eccentricity;   

    fprintf('Area: %.1f px\n', ballArea);
    fprintf('Centroid: (%.2f, %.2f)\n', ballCentroid(1), ballCentroid(2));
    fprintf('Bounding Box: [x=%.2f, y=%.2f, w=%.2f, h=%.2f]\n', ...
        ballBBox(1), ballBBox(2), ballBBox(3), ballBBox(4));
    fprintf('Major axis len: %.2f px\n', ballMajor);
    fprintf('Minor axis len: %.2f px\n', ballMinor);
    fprintf('Perimeter: %.2f px\n', ballPerim);
    fprintf('Eccentricity: %.3f (0 = circle)\n', ballEcc);

    fprintf('Extrema points (x, y):\n');
    for i = 1:size(ballExtrema,1)
        fprintf('  (%.2f, %.2f)\n', ballExtrema(i,1), ballExtrema(i,2));
    end
else
    fprintf('No ball');
end

if any(batMask(:))

    batProps = regionprops(batMask, ...
        'Area', ...
        'Centroid', ...
        'BoundingBox', ...
        'MajorAxisLength', ...
        'MinorAxisLength', ...
        'Perimeter', ...
        'Extrema', ...
        'Orientation');  

    batArea     = batProps.Area;
    batCentroid = batProps.Centroid;
    batBBox     = batProps.BoundingBox;
    batMajor    = batProps.MajorAxisLength;
    batMinor    = batProps.MinorAxisLength;
    batPerim    = batProps.Perimeter;
    batExtrema  = batProps.Extrema;
    batOrient   = batProps.Orientation;     

    fprintf('Area: %.1f px\n', batArea);
    fprintf('Centroid: (%.2f, %.2f)\n', batCentroid(1), batCentroid(2));
    fprintf('Bounding Box: [x=%.2f, y=%.2f, w=%.2f, h=%.2f]\n', ...
        batBBox(1), batBBox(2), batBBox(3), batBBox(4));
    fprintf('Major axis len: %.2f px\n', batMajor);
    fprintf('Minor axis len: %.2f px\n', batMinor);
    fprintf('Perimeter: %.2f px\n', batPerim);
    fprintf('Orientation: %.2f degrees\n', batOrient);

    fprintf('Extrema points (x, y):\n');
    for i = 1:size(batExtrema,1)
        fprintf('  (%.2f, %.2f)\n', batExtrema(i,1), batExtrema(i,2));
    end
else
    fprintf('No bat');
end

figure; imshow(photo); title('Detected Ball (Green) and Bat (Red)');


figure; imshow(photo); title('Detected Ball (Green) and Bat (Red)'); hold on;


ballBounds = bwboundaries(ballMask);
for k = 1:length(ballBounds)
    plot(ballBounds{k}(:,2), ballBounds{k}(:,1), 'g', 'LineWidth', 2);
end


batBounds = bwboundaries(batMask);
for k = 1:length(batBounds)
    plot(batBounds{k}(:,2), batBounds{k}(:,1), 'r', 'LineWidth', 2);
end

hold off;

figure; imshow(photo); title('Original');
figure; imshow(grayPhoto); title('GrayScale');
figure; imshow(photoDouble); title('PhotoDouble');
figure; imshow(grayDouble); title('Gray double');
figure; imshow(grayFiltered); title('Median Filtered photo');
binaryRGB = zeros([size(grayPhoto) 3]);  

binaryRGB(:,:,1) = batMask;   
binaryRGB(:,:,2) = ballMask;  

figure; imshow(binaryRGB); title('Ball (Green) and Bat (Red) on Binary Image');
