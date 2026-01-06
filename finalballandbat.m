clear; clc; close all;
fileNmae = 'TennisSet2.zip';
zipFileName = 'SET2';
unzip(fileNmae, zipFileName);
frames = 21:88;
temp = numel(frames);


ballCentroids = nan(temp,2);   
batCentroids  = nan(temp,2);  
ballSpeed     = nan(temp,1);  
maxBallJump = 80;   
maxBatJump  = 140;  

minBallArea = 10;
maxBallArea = 400;
maxBallEcc  = 0.95;
minBallCirc = 0.35;  


for ii = 1:N
    f = frames(ii);
    fname = sprintf('stennis.%d.ppm', f);
    fullpath = fullfile(zipFileName, fname);

    if ~isfile(fullpath)
        fprintf('Missing frame: %s', fname);
        continue;
    end

    photo = imread(fullpath);
    grayPhoto   = rgb2gray(photo);
    photoDouble = im2double(photo);
    grayFiltered = medfilt2(grayPhoto, [3 3]);

    whiteDet = (photoDouble(:,:,1) > 0.8) & (photoDouble(:,:,2) > 0.8) & (photoDouble(:,:,3) > 0.8);

    lowDetailedPhoto = imgaussfilt(grayFiltered, 5);
    detailsContrasts = (grayFiltered - lowDetailedPhoto) > 0.1;

    ballMaskR = whiteDet & detailsContrasts;
    ballMaskR = imopen(ballMaskR, strel('disk',2));
    ballMaskR = imclose(ballMaskR, strel('disk',2));
    ballMaskR = bwareaopen(ballMaskR, 10);

    [ball, ~] = bwlabel(ballMaskR);
    ballSate = regionprops(ball, 'Area','Centroid','Eccentricity','Perimeter');

    validBall = [];
    for k = 1:numel(ballSate)
        A = ballSate(k).Area;
        P = ballSate(k).Perimeter;
        if P <= 0, continue; end
        circ = 4*pi*A/(P^2);

        if A >= minBallArea && A <= maxBallArea && ...
           ballSate(k).Eccentricity <= maxBallEcc && ...
           circ >= minBallCirc
            validBall(end+1) = k; 
        end
    end

    ballMask = false(size(grayPhoto));

    if ~isempty(validBall)
        prevBall = [NaN NaN];
        for jj = ii-1:-1:1
            if all(~isnan(ballCentroids(jj,:)))
                prevBall = ballCentroids(jj,:);
                break;
            end
        end

        if all(~isnan(prevBall))
            C = vertcat(ballSate(validBall).Centroid);
            d = hypot(C(:,1)-prevBall(1), C(:,2)-prevBall(2));
            [dmin, m] = min(d);

            if dmin <= maxBallJump
                idx = validBall(m);
                ballMask = (ball == idx);
                fprintf('Frame %d: Ball detected (gated). d=%.1f <= %.1f\n', f, dmin, maxBallJump);
            else
                fprintf('Frame %d: Ball NOT found (too far). dmin=%.1f > %.1f\n', f, dmin, maxBallJump);
            end
        else
            areas = [ballSate(validBall).Area];
            [~, m] = max(areas);
            idx = validBall(m);
            ballMask = (ball == idx);
            fprintf('Frame %d: Ball detected (no prior track).\n', f);
        end
    else
        fprintf('Frame %d: No ball candidates after filtering.\n', f);
    end

    redMask   = (photoDouble(:,:,1) > 0.6) & (photoDouble(:,:,2) < 0.4) & (photoDouble(:,:,3) < 0.4);
    blackMask = (photoDouble(:,:,1) < 0.35) & (photoDouble(:,:,2) < 0.35) & (photoDouble(:,:,3) < 0.35);

    batMaskR = redMask | blackMask;
    batMaskR = imopen(batMaskR, strel('rectangle',[3 3]));
    batMaskR = imclose(batMaskR, strel('rectangle',[3 3]));
    batMaskR = bwareaopen(batMaskR, 50);

    [bat, batNums] = bwlabel(batMaskR);
    batStats = regionprops(bat, 'Area', 'MajorAxisLength', 'MinorAxisLength');

    batCandidates = false(size(batMaskR));
    minBatArea = 200; maxBatArea = 2500;
    minBatElongation = 3;

    for k = 1:batNums
        ratio = batStats(k).MajorAxisLength / max(batStats(k).MinorAxisLength, eps);
        if batStats(k).Area >= minBatArea && batStats(k).Area <= maxBatArea && ratio >= minBatElongation
            batCandidates(bat == k) = true;
        end
    end

    [final, ~] = bwlabel(batCandidates);
    propF = regionprops(final, 'Area', 'Centroid');

    batMask = false(size(batMaskR));

    if ~isempty(propF)
        prevBat = [NaN NaN];
        for jj = ii-1:-1:1
            if all(~isnan(batCentroids(jj,:)))
                prevBat = batCentroids(jj,:);
                break;
            end
        end

        if all(~isnan(prevBat))
            C = vertcat(propF.Centroid);
            d = hypot(C(:,1)-prevBat(1), C(:,2)-prevBat(2));
            [dmin, idx] = min(d);

            if dmin <= maxBatJump
                batMask = (final == idx);
            else
                fprintf('Frame %d: Bat NOT found (too far). dmin=%.1f > %.1f\n', f, dmin, maxBatJump);
            end
        else
            [~, idx] = max([propF.Area]);
            batMask = (final == idx);
        end
    else
        fprintf('Frame %d: No bat candidates.\n', f);
    end

    if any(ballMask(:))
        ballProps = regionprops(ballMask, 'Centroid');
        ballCentroid = ballProps.Centroid;
    else
        ballCentroid = [NaN NaN];
    end
    if any(batMask(:))
        batProps = regionprops(batMask, 'Centroid');
        batCentroid = batProps.Centroid;
    else
        batCentroid = [NaN NaN];
    end

    ballCentroids(ii,:) = ballCentroid;
    batCentroids(ii,:)  = batCentroid;

    if ii >= 2
        if all(~isnan(ballCentroids(ii,:))) && all(~isnan(ballCentroids(ii-1,:)))
            dx = ballCentroids(ii,1) - ballCentroids(ii-1,1);
            dy = ballCentroids(ii,2) - ballCentroids(ii-1,2);
            ballSpeed(ii) = hypot(dx, dy);
        end

        if all(~isnan(batCentroids(ii,:))) && all(~isnan(batCentroids(ii-1,:)))
            dx = batCentroids(ii,1) - batCentroids(ii-1,1);
            dy = batCentroids(ii,2) - batCentroids(ii-1,2);
            batSpeed(ii) = hypot(dx, dy);
        end
    end


    if all(~isnan(ballCentroid))
        fprintf('Ball centroid: (%.2f, %.2f)\n', ballCentroid(1), ballCentroid(2));
    else
        fprintf('Ball centroid not found\n');
    end

    if all(~isnan(batCentroid))
        fprintf('Bat  centroid: (%.2f, %.2f)\n', batCentroid(1), batCentroid(2));
    else
        fprintf('Bat  centroid nit found\n');
    end

    if ii >= 2
        if ~isnan(ballSpeed(ii))
            fprintf('Ball speed (frame %d -> %d): %.3f px/frame\n', frames(ii-1), frames(ii), ballSpeed(ii));
        else
            fprintf('Ball speed (frame %d -> %d): N/A\n', frames(ii-1), frames(ii));
        end

        if ~isnan(batSpeed(ii))
            fprintf('Bat  speed (frame %d -> %d): %.3f px/frame\n', frames(ii-1), frames(ii), batSpeed(ii));
        else
            fprintf('Bat  speed (frame %d -> %d): N/A\n', frames(ii-1), frames(ii));
        end
    else
        fprintf('Speed: N/A (first frame)\n');
    end

    figure(100); clf;
    imshow(photo);
    title(sprintf('Frame %d: Detected Ball (Green) and Bat (Red)', f));
    hold on;

    ballBoundaries = bwboundaries(ballMask);
    for k = 1:length(ballBoundaries)
        plot(ballBoundaries{k}(:,2), ballBoundaries{k}(:,1), 'g', 'LineWidth', 2);
    end

    batBoundaries = bwboundaries(batMask);
    for k = 1:length(batBoundaries)
        plot(batBoundaries{k}(:,2), batBoundaries{k}(:,1), 'r', 'LineWidth', 2);
    end

    hold off;
    drawnow;
end


fprintf('Frame | BallCx  BallCy | BatCx   BatCy  | BallSpd | BatSpd\n');
for ii = 1:N
    fprintf('%5d | %6.2f %6.2f | %6.2f %6.2f | %7.3f | %7.3f\n', ...
        frames(ii), ...
        ballCentroids(ii,1), ballCentroids(ii,2), ...
        batCentroids(ii,1),  batCentroids(ii,2), ...
        ballSpeed(ii), batSpeed(ii));
end

figure;
plot(frames, ballSpeed, '-o', 'LineWidth', 1.5); hold on;
plot(frames, batSpeed,  '-o', 'LineWidth', 1.5);
grid on;
xlabel('Frame number');
ylabel('Speed (pixels/frame)');
title('2D Velocity Magnitude (Ball and Bat)');
legend('Ball', 'Bat', 'Location', 'best');


binaryRGB = zeros([size(grayPhoto) 3]);
binaryRGB(:,:,1) = batMask;   
binaryRGB(:,:,2) = ballMask;  

figure; imshow(binaryRGB);
title('Ball (Green) and Bat (Red) on Binary Image (Last Frame Processed)');
