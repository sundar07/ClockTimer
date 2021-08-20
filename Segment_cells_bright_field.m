function Segment_cells_bright_field()

% Created by Sundar Naganathan - 2016-2017
% This code segments single cells from brightfield images. The code is set for detecting cells that are approx. at the centre of the image, but can be
% adapted to segment cells anywhere in the image

movie_path = cd;

% Make a new folder called 'Seg_images' if it does not already exist
if exist(strcat(movie_path,'/Seg_images'),'dir') ~= 7
   mkdir(strcat(movie_path,'/Seg_images'))
end

% Get tif file information
DirOutput = dir(strcat(movie_path,'/*.tif'));
FileNames = {DirOutput.name};
FileInfo = imfinfo(FileNames{1});

% Pre-initialize the variable
I1 = nan(FileInfo(1).Height,FileInfo(1).Width,numel(FileInfo));

% Read all frames
for Frame = 1:numel(FileInfo)
    I1(:,:,Frame) = im2double(imread(strcat(movie_path,'/',FileNames{1}),Frame)); % DIC
end

Segment.mfile_used = 'Segment_cells_bright_field';

% Get center position of the image
Xcenter = size(I1,2)/2;
Ycenter = size(I1,1)/2;

Total_num_frames = size(I1,3);

% Pre-initialize more variables
Segment.Mask = nan(size(I1));
Segment.Area = nan(Total_num_frames,1);
Segment.Centroid = nan(Total_num_frames,2);

for frame = 1:Total_num_frames

    Iframe1 = I1(:,:,frame); % choose the first frame

    Iad = adapthisteq(Iframe1); % Enhance the contrast of the image

    % Filter the image by considering local pixel values. I have chosen a really small value for 'DegreeOfSmoothing' so as to ensure there 
    % is no smoothing around edges of cells
    Ifilt = imguidedfilter(Iad,'NeighborhoodSize',[3 3],'DegreeOfSmoothing',0.001);

    SE = strel('disk',2); % Define a structuring element
    Ierode = imerode(Ifilt,SE); % Erode the image using the structural element
    Idilate = imdilate(Ifilt,SE); % Dilate the image using the structural element
    Igrad = Idilate - Ierode; % Generate gradient by subtracting eroded image from the dilated image

    % Threshold and fill holes in image
    Ithresh = imbinarize(Igrad,graythresh(Igrad));

    Iareaopen = bwareaopen(Ithresh,300); % I give a large value for connected components as I dont want anything apart from a segmented cell
    Ilabel = bwlabel(Iareaopen); % Label all detected components
    R = regionprops(Ilabel); % Get statistics of all labeled components
    
    % If the number of objects are lesser than an arbitrary number (10 seems to be a good number for our images), then a further round of erosion is not needed
    % If its greater than this arbitrary number, then subjecting the images to an erosion gives better results
    if numel(R) <= 10
        Ifill = imfill(Iareaopen,'holes'); % Fill holes
        Iclearborder = imclearborder(Ifill); % clear any segments attached to the border
        Idilate = imdilate(Iclearborder,strel('disk',3)); % dilate the image a little bit so that you can use contraction-based algorithm in activecontour
        Iareaopen2 = bwareaopen(Idilate,300); % remove smaller segmented objects
    else
        Ierode = imerode(Iareaopen,strel('disk',4)); % Do an erosion 
        Iclearborder = imclearborder(Ierode); % clear any segments attached to the border
        Idilate = imdilate(Iclearborder,strel('disk',8)); % dilate the image a little bit so that you can use contraction-based algorithm in activecontour
        Ifill2 = imfill(Idilate,'holes'); % Fill holes
        Iareaopen2 = bwareaopen(Ifill2,300); % remove smaller segmented objects
    end

    % Get area and centroid information about segmented regions
    Ilabel2 = bwlabel(Iareaopen2);
    R2 = regionprops(Ilabel2,'Area','Centroid');
    area = [R2.Area];
    centroids = cat(1,R2.Centroid);

    % If nothing was detected, then proceed to next frame
    if ~isempty(centroids)

        % Find distance of segmented regions from image center
        Dist_from_center = sqrt((centroids(:,1)-Xcenter).^2 + (centroids(:,2)-Ycenter).^2);
        [~,dist_idx] = sort(Dist_from_center,'ascend');

        [~,area_idx] = sort(area,'descend');  
        
        % If the number of segmented regions are more than 1 but less than 5, use the same idx as area_idx
        % elseif the number of segmented regions are more than 5, find common regions by considering just the first five regions from both area_idx and dist_idx
        % If number of segmented regions is equal to 1, then idx_cell is 1
        if numel(area_idx) > 1 && numel(area_idx) <= 5
            idx_cell = area_idx;
        elseif numel(area_idx) > 5
            idx_cell = intersect(area_idx(1:5),dist_idx(1:5),'stable');   
        else
            idx_cell = 1;
        end
        
        % Just in case, no common regions are found, increase the number of regions to be compared to 10
        % Usually at least one common region will be found. If no common regions are still found, then segmentation did not work properly and in the next step you proceed to skip this frame
        if isempty(idx_cell)
            if numel(area_idx) > 1 && numel(area_idx) <= 10
                idx_cell = area_idx;
            elseif numel(area_idx) > 10
                idx_cell = intersect(area_idx(1:10),dist_idx(1:10),'stable');   
            else
                idx_cell = 1;
            end
        end

        % If more than one region is identified, choose the one that has the highest standard deviation of pixel intensities 
        % and from second frame onwards, do an additional check with distance from segmented cell from previous frame. Usually cells do not move much, 
        % so the one closest to the previous frame should be the cell of interest
        % If no region is identified, then break the code and move to the next frame
        if numel(idx_cell) == 1
            Icell = ismember(Ilabel2,idx_cell);
            Cell_interest = idx_cell;
        elseif isempty(idx_cell)
            % Alternatively you can add a little piece of code here prompting the user to select the region of interest manually, which can then be subjected to activecontour algorithm
            disp(['Cell detection for frame ' num2str(frame) 'did not work properly'])
            continue;
        elseif numel(idx_cell) > 1        

            STD_regions = nan(numel(idx_cell),1); % Pre-initialize the variable
            % Get standard deviation of pixel intensities of each segmented region
            for ii = 1:numel(idx_cell)
                Isel = Ilabel2 == idx_cell(ii);
                Imult = Isel .* Iframe1;
                STD_regions(ii,1) = nanstd(Imult(:));
            end
            % the one with the maximum standard deviation is usually the cell
            [~,Max_STD_regions_idx] = max(STD_regions);
            Cell_interest = idx_cell(Max_STD_regions_idx);
            Icell = ismember(Ilabel2,Cell_interest);
            
            % If frame is greater than 1, then proceed to check with centroid of cell from previous frame
            if frame > 1
                Centroid_seg_regions = centroids(idx_cell,1:2); % Get all centroids
                % Find distance between cell from previous frame and all segmented regions from current frame
                Dist_from_previous_cell = sqrt((Centroid_seg_regions(:,1)-Segment.Centroid(frame-1,1)).^2 + (Centroid_seg_regions(:,2)-Segment.Centroid(frame-1,2)).^2);
                
                % Choose the region with maximum std and get the distance of this region from cell from previous frame
                Dist_cell_interest_from_previous = Dist_from_previous_cell(idx_cell == Cell_interest,1);
                
                % If distance is less than 50, then the one with max std is the cell of interest, so Icell from line 139 is the right one
                % If distance is more than 50, then find the closest one from previous cell and use that instead of Icell from line 139
                if Dist_cell_interest_from_previous > 50
                    [~,Least_dist_idx] = min(Dist_from_previous_cell);
                    Cell_interest = idx_cell(Least_dist_idx);
                    Icell = ismember(Ilabel2,Cell_interest);
                end
            end
        end

        % Apply active contour algorithm to segment cell
        % A positive ContractionBias makes sure that the segmented region gets constricted while applying the algorithm
        Icontour = activecontour(Igrad,Icell,300,'Chan-Vese','SmoothFactor',1,'ContractionBias',0.1);

        % Sometimes, this algorithm might segment the middle region into 2 or more regions
        % Therefore get number of segmented regions and choose the biggest region, which will be the cell
        Ilabel2 = bwlabel(Icontour);
        R2 = regionprops(Ilabel2,'Area');
        area2 = [R2.Area];
        [~,idx3] = max(area2);

        Icell2 = ismember(Ilabel2,idx3);

        R3 = regionprops(Icell2,'Area','Centroid');
        % If even after all these steps, more than 1 region is identified, then segmentation did not work properly and move on to next frame
        % Alternatively you can add a little piece of code here prompting the user to select the region of interest
        if numel(R3) == 1
            Segment.Area(frame,1) = [R3.Area];
            Segment.Centroid(frame,1:2) = cat(1,R3.Centroid);
        else
            disp(['Cell detection for frame ' num2str(frame) 'did not work properly'])
            continue;
        end

        Mask = nan(size(Iframe1));
        Mask(Icell2 == 1) = 1;
        Segment.Mask(:,:,frame) = Mask;

        % Get regions of interest around the cell for determining background
        Bg1_2_x = [Segment.Centroid(frame,1)-80,Segment.Centroid(frame,1)-80,Segment.Centroid(frame,1)-50,Segment.Centroid(frame,1)-50,Segment.Centroid(frame,1)-80];
        Bg1_4_y = [Segment.Centroid(frame,2)-80,Segment.Centroid(frame,2)-50,Segment.Centroid(frame,2)-50,Segment.Centroid(frame,2)-80,Segment.Centroid(frame,2)-80];
        Bg2_3_y = [Segment.Centroid(frame,2)+50,Segment.Centroid(frame,2)+80,Segment.Centroid(frame,2)+80,Segment.Centroid(frame,2)+50,Segment.Centroid(frame,2)+50];
        Bg3_4_x = [Segment.Centroid(frame,1)+50,Segment.Centroid(frame,1)+50,Segment.Centroid(frame,1)+80,Segment.Centroid(frame,1)+80,Segment.Centroid(frame,1)+50];

        Segment.Bg1_mask(:,:,frame) = roipoly(Mask,Bg1_2_x,Bg1_4_y);
        Segment.Bg2_mask(:,:,frame) = roipoly(Mask,Bg1_2_x,Bg2_3_y);
        Segment.Bg3_mask(:,:,frame) = roipoly(Mask,Bg3_4_x,Bg2_3_y);
        Segment.Bg4_mask(:,:,frame) = roipoly(Mask,Bg3_4_x,Bg1_4_y);

        % Get cell and background masks in a single image, which can be used for display purposes
        Cell_and_background = zeros(size(Iframe1));
        Cell_and_background(Icell2 == 1) = 1;
        Cell_and_background(Segment.Bg1_mask(:,:,frame) == 1) = 1;
        Cell_and_background(Segment.Bg2_mask(:,:,frame) == 1) = 1;
        Cell_and_background(Segment.Bg3_mask(:,:,frame) == 1) = 1;
        Cell_and_background(Segment.Bg4_mask(:,:,frame) == 1) = 1;

        % Get perimeter of cell and backgrounds
        Iperim = bwperim(Cell_and_background);
        % Overlay segmented image perimeter on top of original image
        Ioverlay = imoverlay(Iframe1,Iperim,[1 0 0]);
        % Display segmented image and save it
        figure(1);clf;imshow(Ioverlay)
        saveas(gcf,strcat(movie_path,'/Seg_images/Image_bg_sub_',num2str(frame)),'tif')


    end % if ~isempty(centroids)
end % for frame = 1:End_frame_analysis

save(strcat(movie_path,'/Segment_bg_sub.mat'),'Segment')