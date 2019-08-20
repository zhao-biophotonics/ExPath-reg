%Rigid registration of pre-expansion and post-expansion images using SIFT and RANSAC.
%Created by Yongxin based on the demo code in vlfeat.org 20151208
%Added downsampling to speed up the search, by Yongxin Zhao 20151210
%Generalize the z projection matching, by Yongxin zhao 20160404
%Speed up the process by skipping adjacent z plane if the number of matching keypoints is lower than the user defined value, by Yongxin Zhao 20180411
%Code clean-up by Yongxin Zhao 20190820
%The code uses vlfeat for SIFT and RANSAC computation.
% @misc{vedaldi08vlfeat,
%  Author = {A. Vedaldi and B. Fulkerson},
%  Title = {{VLFeat}: An Open and Portable Library
%           of Computer Vision Algorithms},
%  Year  = {2008},
%  Howpublished = {\url{http://www.vlfeat.org/}}
% }
%The input are a pair of tiff images. The output is a folder with
%registerred images and some statistics, including expansion factor.
%Autoright

vl_setup
close all;
clearvars -except path_name;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustable parameters
tileSize = 27000; %Decrease if out of memory
peakThresh = 0; %SIFT peak threshold
edgeThresh = 6; %SIFT edge threshold
nRansacTrials = 40000; %Increase for more reliable matching
nPtsFit = 5; %For each RANSAC trial
nKeypointsThresh = 17; %Minimum number of keypoint matches
radius = 30; %Maximum tolerated RANSAC distance
ChannelNumber = 4; %The number of channels for the pre-expansion image (map)
ChannelNumber2 = 4; %The number of channels for the post-expansion image (query)
z_pre = 1; %The number of planes for maximum intensity projection for pre-expansion images. 1 means no z projection.
z_projection =5; %The number of planes for maxium intensity projection in the query images. 1 means no z projection.
downsampling_pre = 2; %Downsampling factor to speed up the search
downsampling_post = 10; %Downsampling factor to speed up the search
map_pixel = 0.163; % Pixel size of the map image, unit: um. 
query_pixel = 0.163; %Pixel size of the query image, unit: um.
keypoint_threshold = 7; %User-defined threshold for skipping adjacent z plane.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loading images.
path_name='G:\Dropbox (MIT)\Imaging\ExM\20170531 Breast lesions\'; 
mapName='abc.tif';
% if exist('mapName')
%     if isequal(mapName,0)
%         [mapName, path_name] = uigetfile({'*.tif';'*.png';'*.*'},'Find the map image');
%     else
        [mapName, path_name] = uigetfile({'*.tif';'*.png';'*.*'},'Find the map image',[path_name,mapName]);
%     end
%     else
%     [mapName, path_name] = uigetfile({'*.tif';'*.png';'*.*'},'Find the map image');
% end
[queryName, path_name1] = uigetfile({'*.tif';'*.png';'*.*'},'Find the query image',[path_name,mapName]);
%%
map_data = TIFFStack([path_name,mapName]);
[~,~,map_size] = size(map_data);

%Loading query data
query_data = TIFFStack([path_name1,queryName]);
[query_x,query_y,query_size] = size(query_data);

disp(['Loaded map image ',mapName]);
disp(['Loaded query image ',queryName]);

%Initiate parameters
keypoints = zeros (map_size, query_size);

for i=1:ChannelNumber
    for ii=i:ChannelNumber2:(query_size+ChannelNumber2*(1-z_projection))
    ii_seq=[];%empty the variable ii_seq for every cycle.
        for iiz=1:z_projection %To generate mean intensity projection
            iie = ii + (iiz-1)*ChannelNumber2;
            ii_seq=[ii_seq iie];
        end
        for iii=i:ChannelNumber:(map_size+ChannelNumber*(1-z_pre))
            if iii==i
            else
             if (keypoints(iii-ChannelNumber,ii)>0 && keypoints(iii-ChannelNumber,ii)< keypoint_threshold)
                continue
             end
            end
            iii_seq=[];
            for iiiz=1:z_pre %To generate mean intensity projection
            iiie = iii + (iiiz-1)*ChannelNumber;
            iii_seq=[iii_seq iiie];
            end
            mapOrig = map_data(:,:,iii_seq);
            mapOrig = uint16(mean(mapOrig,3));
            map = imadjust(im2single(imresize(mapOrig,1/downsampling_pre)));
            queryOrig = query_data(:,:,ii_seq);
            queryOrig = uint16(mean(queryOrig,3));
            query = im2single(imadjust(imresize(queryOrig,1/downsampling_post)));
            
            %Subtract background
            map_avg = mean(map(:));
            map_std = std(map(:));
            map_threshold = single (map_avg - 0.5*map_std);
%             map_mask = map >=map_threshold;
            map = map -map_threshold;
            map (map(:)<0)=0;
            
            query_avg = mean(query(:));
            query_std = std(query(:));
            query_threshold = single (query_avg - 0.5*query_std);
%             query_mask = query >=query_threshold;
%             query = query.*query_mask;
            query = query -query_threshold;
            query (query(:)<0)=0;
            
            fprintf('Checking: Channel %d, map z number: %d, query z number %d...',...
    i, fix((iii-1)/ChannelNumber)+1,fix((ii-1)/ChannelNumber2)+1);
            [~,~,keypoints(iii,ii)]=roifinder2(map,query);
            %Add conditions to speed up the process, if too few keypoints,
            %jump one or two z slices.
                                  
           
        end
    end
end

[max_keypoint,index_keypoint]=max(keypoints(:));
[map_best, query_best] = ind2sub(size(keypoints),index_keypoint);
bestChannel=mod(map_best, ChannelNumber);
if bestChannel==0
    bestChannel=ChannelNumber;
end

fprintf('Best match: Channel %d, map z number: %d, query z number %d, with %d keypoints matched\n\n',...
    bestChannel, fix((map_best-1)/ChannelNumber)+1,fix((query_best-1)/ChannelNumber)+1,max_keypoint);
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%crop the part that matches the query images. 
defaultChannel = bestChannel; %The default channel for matching
%To determine the lower and upper range of the best z plane for the map
%image.
lowrange_map = fix((map_best-1)/ChannelNumber)*ChannelNumber+1;
uprange_map = (fix((map_best-1)/ChannelNumber)+1)*ChannelNumber;
%This part is for the query image.
lowrange_query = fix((query_best-1)/ChannelNumber2)*ChannelNumber2+1;
uprange_query = (fix((query_best-1)/ChannelNumber2)+1)*ChannelNumber2;
%Extract the best image plane for both stacks.
%map_best_data = map_data (:,:,lowrange_map:uprange_map);
%mapOrig = map_best_data(:,:,defaultChannel);

for i=1:ChannelNumber
    iii_seq=[];%empty the variable ii_seq for every cycle.
    ii=lowrange_map+i-1;
    for iiiz=1:z_pre %To generate mean intensity projection
            iiie = ii + (iiiz-1)*ChannelNumber;
            iii_seq=[iii_seq iiie];
    end
    mapOrig = map_data(:,:,iii_seq);
    mapOrig = uint16(mean(mapOrig,3));
        if exist('mapOrig_mean')
            mapOrig_mean = cat(3,mapOrig_mean,mapOrig);
        else
            
            mapOrig_mean = mapOrig;
        end
    
end
map = im2single(imadjust(mapOrig_mean(:,:,defaultChannel)));
            %Subtract background
            map_avg = mean(map(:));
            map_std = std(map(:));
            map_threshold = single (map_avg - 0.5*map_std);
%             map_mask = map >=map_threshold;
            map = map -map_threshold;
            map (map(:)<0)=0;
%Get the best z projection for the post-expansion images
%%
for i=1:ChannelNumber2
    ii_seq=[];%empty the variable ii_seq for every cycle.
    ii=lowrange_query+i-1;
    for iiz=1:z_projection %To generate mean intensity projection
            iie = ii + (iiz-1)*ChannelNumber2;
            ii_seq=[ii_seq iie];
    end
    queryOrig = query_data(:,:,ii_seq);
    queryOrig = uint16(mean(queryOrig,3));
        if exist('queryOrig_mean')
            queryOrig_mean = cat(3,queryOrig_mean,queryOrig);
        else
            
            queryOrig_mean = queryOrig;
        end
    
end

query = im2single(imadjust((imresize(queryOrig_mean(:,:,defaultChannel),1/downsampling_post*downsampling_pre))));


            query_avg = mean(query(:));
            query_std = std(query(:));
            query_threshold = single (query_avg - 0.5*query_std);
            query = query -query_threshold;
            query (query(:)<0)=0;



% Generate database of tiles
[m,n] = size(map);
N = ceil(n/tileSize);
M = ceil(m/tileSize);
nTiles = M*N;
database = cell(1,nTiles);
for i = 1:M
    for j = 1:N
        img = im2single(map((i-1)*floor(m/M)+1:i*min(floor(m/M),m),...
            (j-1)*floor(n/N)+1:j*min(floor(n/N),n)));
        database{(i-1)*N+j} = img/max(img(:)); %Normalize intensity
    end
end

% % Check if SIFT keypoints already exist in workspace
% isPrecomputed = exist('f1s','var') && exist('d1s','var') ...
%     && exist('database','var') && iscell(f1s) && iscell(d1s) ...
%     && iscell(database) && (numel(f1s) == numel(database)) && ...
%     (numel(d1s) == numel(database)) && isequal(database,db_old);

% Compute SIFT keypoints for map
% if isPrecomputed
%     disp('Existing SIFT keypoints found in workspace.');
% else
    f1s = cell(1,nTiles);
    d1s = cell(1,nTiles);
    disp('Computing SIFT keypoints...');
    reverseStr = '';
    for i = 1:nTiles
        [f1s{i},d1s{i}] = vl_sift(database{i},...
            'PeakThresh',peakThresh,'EdgeThresh',edgeThresh);
        msg = sprintf('Completed tile %d of %d\n',i,nTiles);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
    end
% end

% SIFT keypoints for query image
[f2,d2] = vl_sift(query,'PeakThresh',peakThresh,...
    'EdgeThresh',edgeThresh);

% Find matching image
nMatches = zeros(1,nTiles);
matchResult = cell(1,nTiles);
disp('Matching with isometry model...');
reverseStr = '';
for i = 1:nTiles
    % SIFT matches
    tile = database{i};
    f1 = f1s{i};
    d1 = d1s{i};
    [matches, ~] = vl_ubcmatch(d1,d2); % Distance ratio test
    
    % Remove many-to-one matches
    [uniqueRow2, IA, ~] = unique(matches(2,:));
    uniqueRow1 = matches(1,IA);
    matches = [uniqueRow1; uniqueRow2];
    numMatches = size(matches,2);
    
    X1 = f1(1:2,matches(1,:)); X1(3,:) = 1;
    X2 = f2(1:2,matches(2,:)); X2(3,:) = 1;
    
    % RANSAC with geometric model
    H = cell(1,nRansacTrials);
    isMatch = cell(1,nRansacTrials);
    score = zeros(1,nRansacTrials);
    matched_mse = zeros(1,nRansacTrials);
    for j = 1:nRansacTrials
        % Estimate isometry
        subset = vl_colsubset(1:numMatches,nPtsFit);
        X = X1(1:2,subset);
        Y = X2(1:2,subset);
        [Q,s,t] = fit_isometry(X,Y);
        H{j} = [cat(2,s*Q,t); 0 0 1];
        % Score isometry
        X2_ = H{j}*X1;
        du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:);
        dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:);
        isMatch{j} = (du.^2 + dv.^2 < radius^2);
        score(j) = sum(isMatch{j});
        matched_mse(j) = sum((du.^2 + dv.^2).*isMatch{j});
    end
    
    % Find best mapping for current tile
    [~,best] = find(score == max(score));
    [~,idx] = min(matched_mse(best));
    H = H{best(idx)};
    isMatch = isMatch{best(idx)};
    msg = sprintf('Completed tile %d of %d\n',i,nTiles);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'),1,length(msg));
    
    % Plot
    figure(1); subplot(M,N,i); hold on;
    dh1 = max(size(query,1)-size(tile,1),0);
    dh2 = max(size(tile,1)-size(query,1),0);
    imshow([padarray(tile,dh1,'post') padarray(query,dh2,'post')],...
        'InitialMagnification','fit');
    title(sprintf('Tile %d of %d\n%d keypoints matched',...
        i,nTiles,sum(isMatch)));
    o = size(tile,2);
    hold on;
    plot([f1(1,matches(1,isMatch))+1;f2(1,matches(2,isMatch))+o+1], ...
        [f1(2,matches(1,isMatch))+1;f2(2,matches(2,isMatch))+1],...
        ['y' '-'],'LineWidth',1);
    plot(f1(1,matches(1,isMatch))+1,f1(2,matches(1,isMatch))+1,...
        'yo','MarkerSize',5,'MarkerFaceColor','y');
    plot(f2(1,matches(2,isMatch))+o+1,f2(2,matches(2,isMatch))+1,...
        'ko','MarkerSize',5,'MarkerFaceColor','y');
    hold off;
    newfolder = [path_name1, queryName(1:end-4),'_ds',num2str(downsampling_post)];
    mkdir (newfolder);
    saveas (gca, [newfolder,'\','key points.fig']);
    % Update
    matchResult{i}.ratio_test = size(matches,2);
    matchResult{i}.ransac = sum(isMatch);
    matchResult{i}.mse = matched_mse(best(idx));
    matchResult{i}.model = H;
    matchResult{i}.matches = matches(:,isMatch);
    matchResult{i}.f1 = f1;
    matchResult{i}.f2 = f2;
    matchResult{i}.loc1 = f1(1:2,matches(1,isMatch));
    matchResult{i}.loc2 = f2(1:2,matches(2,isMatch));
    
    nMatches(i) = sum(isMatch);
    if nMatches(i) >= nKeypointsThresh
        disp('Match found!');
        break;
    end
end

% Find best-matched tile
[~,idx] = max(nMatches);
matchResult = matchResult{idx};
fprintf('Best match: Tile %d of %d with %d keypoints matched\n\n',...
    idx,nTiles,max(nMatches));

% Calculate image borders using homography
H = inv(matchResult.model);
[h,w] = size(query);
x = [(1:w)'; (1:w)'; ones(h-2,1); w*ones(h-2,1)];
y = [ones(w,1); h*ones(w,1); (2:h-1)'; (2:h-1)'];
z = H(3,1)*x + H(3,2)*y + H(3,3);
x_ = (H(1,1)*x + H(1,2)*y + H(1,3))./z;
y_ = (H(2,1)*x + H(2,2)*y + H(2,3))./z;
loc = [x_,y_]';

% Display
disp(matchResult);
nLoc = size(loc,2);
mask = zeros(size(map));
i = ceil(idx/M);
j = idx - (i-1)*N;
x_offset = (i-1)*floor(m/M);
y_offset = (j-1)*floor(n/M);
y = round(loc(1,:)) + y_offset*ones(1,nLoc);
x = round(loc(2,:)) + x_offset*ones(1,nLoc);
for k = 1:nLoc
    if 1 <= y(k) && y(k) <= n && 1 <= x(k) && x(k) <= m
        mask(x(k),y(k)) = 1;
    end
end
warning('off','Images:initSize:adjustingMag');
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(121); imshow(imadjust(im2double(mapOrig) + ...
    imdilate(mask,strel('disk',8))));
subplot(122); imshow(imtransform(query,maketform('affine',H'),'XYScale',1));
warning('on','Images:initSize:adjustingMag');

% Crop all the channel with the boundary and save them all in a new folder.

mask_backup = mask;
for ii=1:ChannelNumber
    mask = mask_backup; % reset the mask in every cycle
    ly = find(sum(mask) > 0,1,'first');
    ry = find(sum(mask) > 0,1,'last');
    ux = find(sum(mask,2) > 0,1,'first');
    lx = find(sum(mask,2) > 0,1,'last');
    cropped_roi = mapOrig_mean(ux:lx,ly:ry, ii);
    mask = mask(ux:lx,ly:ry);

    H = [cat(2,normc(matchResult.model(1:2,1:2))',[0;0]); 0 0 1];
    tform = maketform('affine',H);
    cropped_roi = imtransform(cropped_roi,tform,'XYScale',1);
    mask = imtransform(mask,tform);
    ly = find(sum(mask) > 0,1,'first');
    ry = find(sum(mask) > 0,1,'last');
    ux = find(sum(mask,2) > 0,1,'first');
    lx = find(sum(mask,2) > 0,1,'last');
    cropped_roi = cropped_roi(ux:lx,ly:ry);
    imwrite(cropped_roi,[newfolder '\','C', num2str(ii), '_cropped_pre_best_', mapName]);
    
end

for ii=1:ChannelNumber2
    imwrite(queryOrig_mean(:,:,ii),[newfolder '\','C', num2str(ii), '_cropped_post_best_', queryName]);
end
%calculate expansion factors, assuming all the imaging parameters are the
%same.
[cropped_x,cropped_y] = size (cropped_roi);
expansion_factor = mean(query_x/cropped_x/2+query_y/cropped_y/2)*query_pixel/map_pixel;
fileID = fopen([newfolder '\note.txt'],'w');
fprintf(fileID, 'path name: %s\r\n',path_name);
fprintf(fileID, 'map file name: %s\r\n',mapName);
fprintf(fileID, 'query file name: %s\r\n',queryName);
fprintf(fileID,'Best match: Channel %d, map z number: %d, query z number %d, with %d keypoints matched\r\n',...
    bestChannel, fix((map_best-1)/ChannelNumber)+1,fix((query_best-1)/ChannelNumber)+1,max_keypoint);
fprintf(fileID, 'map pixel size: %s\r\n',map_pixel);
fprintf(fileID, 'query pixel size: %s\r\n',query_pixel);
fprintf(fileID, 'The expansion factor is %d\n',expansion_factor);
fclose(fileID);
saveas (gca, [newfolder,'\','cropping result.fig']);
save([newfolder '\',  'variable.mat']);
