function [thetas,nosem]=GetWhiskerAngle(vid,angleRectify,whiskerRegion,noseRegion)

nframes = floor(vid.Duration * vid.FrameRate);

%%% some initialization %%%
maxNumFramesPerBatch = 10000; % 10000;
numBatches = ceil(nframes / maxNumFramesPerBatch);
thetas = zeros(1, nframes);
nosem = zeros(2, nframes);

% stride=2;

for batchNum = 1: numBatches
    tic
    %%% read frames %%%
    range = [(batchNum - 1) * maxNumFramesPerBatch + 1, min(batchNum * maxNumFramesPerBatch, nframes)];
    data = zeros(whiskerRegion(4)-whiskerRegion(3)+1,whiskerRegion(2)-whiskerRegion(1)+1, diff(range) + 1);
    for i = 1 : diff(range) + 1
        try 
        tmp = double(rgb2gray(readFrame(vid)));
        catch ME
            if strcmp(ME.identifier,'MATLAB:audiovideo:VideoReader:EndOfFile')
                disp('No more frames available to read from file.')
            end
            continue
        end
%         tmp = imresize(tmp, scale);
        tmp = imrotate(tmp, angleRectify);
        data(:, :, i) = tmp(whiskerRegion(3):whiskerRegion(4),whiskerRegion(1): whiskerRegion(2));
        if mod(i, 100) == 0
            disp(['Done frames ', num2str(i), '/', num2str(nframes), ' - batch ', num2str(batchNum)])
        end
    end
    if sum(sum(sum(data)))==0
        continue
    end
    data = WT_normalize(data);
    toc
    %%% get nose info %%%
    nosec = zeros(2, diff(range) + 1);
    datanose = data(noseRegion(3):noseRegion(4),noseRegion(1): noseRegion(2), :);
    mxnose = max(datanose, [], 3);
    datanose = mxnose - datanose;
    datanose = WT_normalize(datanose);
    sigma=3;
    for i = 1: diff(range) + 1
        tmp = datanose(:, :, i);
        tmp = imgaussfilt(tmp, sigma);
        tmp = tmp > graythresh(tmp);
        tt = regionprops(tmp, 'Centroid');
        while length(tt)>1
         tmp = imgaussfilt(double(tmp), sigma);
         tmp = tmp > graythresh(tmp);
         tt = regionprops(tmp, 'Centroid');
         sigma=sigma+1;
%          figure(gcf);imagesc(tmp)
        end
        nosec(:,i) = tt.Centroid';
    end
    nosem(:,range(1): range(2)) = nosec;
    
    %%% get mask, optional %%%
    mx = max(data, [], 3);
    mxt = imgaussfilt(mx, 10);
    hd = mxt < 0.1;
    hdt = imdilate(hd, strel('disk', 11));
    mask = 1 - hdt;
    [l, n] = bwlabeln(mask);
    for i = 1: n
        tt = l == i;
        s(i) = sum(tt(:));
    end
    [~, i] = max(s);
    mask = l == i;
    
    %%% preprocess %%%
%     ispara = 1;
%     dt1 = WT_anidenoise(data, ispara);
    dt1 = 1 - WT_normalize(data);
    dt2 = WT_bg_remove(dt1, 5);
    datan = zeros(size(dt2));
    datan(6: end - 5, 6: end - 5, :) = dt2(6: end - 5, 6: end - 5, :);
    datan = datan .* mask;
    clear data
    clear dt1
    clear dt2
    
    %%% get whisking angle %%%
    thetac = zeros(1, diff(range) + 1);
    for i = 1: diff(range) + 1
        tmp = datan(:, :, i);
        img = tmp > graythresh(tmp);
        [H, theta, rho] = hough(img);
        P = houghpeaks(H,10,'threshold',ceil(0.3*max(H(:))));
        lines = houghlines(img,theta,rho,P,'FillGap',5,'MinLength',7);
        idx = zeros(1, length(lines));
        for j = 1: length(lines)
            idx(j) = lines(j).theta;
        end
        idn = find(idx < 0);
        idx(idn) = 180 + idx(idn);
        thetac(i) = median(idx);
    end
    
    thetas(range(1): range(2)) = thetac;
    disp(['Done batch # ', num2str(batchNum)])
end