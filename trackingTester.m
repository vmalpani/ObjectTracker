function trackingTester(data_params, tracking_params)
    kernel = computeEpanechnikovKernel(tracking_params.rect(3),tracking_params.rect(4));
    
    % reading images into a stack
    img_stack = cell(length(data_params.frame_ids));
    for i = data_params.frame_ids
        img_stack{i} = im2double(imread(fullfile(data_params.data_dir,...
           data_params.genFname(data_params.frame_ids(i)))));
    end
    
    % all computations for first image outside the for loops 
    % to save the color map to be used by subsequent frames
    img1 = img_stack{1};
    % create output directory
    mkdir (fullfile(data_params.out_dir));
    img1 = drawBox(img1, tracking_params.rect, [0, 0, 255], 3);
    imwrite(img1, fullfile(data_params.out_dir,...
    data_params.genFname(data_params.frame_ids(1))));
    
    % median filtering image to smooth out noise
    smooth_img1(:,:,1) = medfilt2(img1(:,:,1),[tracking_params.filter_size tracking_params.filter_size]);
    smooth_img1(:,:,2) = medfilt2(img1(:,:,2),[tracking_params.filter_size tracking_params.filter_size]);
    smooth_img1(:,:,3) = medfilt2(img1(:,:,3),[tracking_params.filter_size tracking_params.filter_size]);
    
    % crop the object of interest
    bounded_object = smooth_img1(tracking_params.rect(2):tracking_params.rect(2)+tracking_params.rect(4),tracking_params.rect(1):tracking_params.rect(1) + tracking_params.rect(3),:);
    
    figure, imshow(bounded_object);
    % generating index map and histogram
    [ind1,map1] = rgb2ind(bounded_object,tracking_params.bin_n);
    img1_hist = imhist(ind1(:), map1);

    % reading images as per the ids specified
    for i = data_params.frame_ids(1)+1:data_params.frame_ids(end) 
       img = img_stack{i};
       % median filtering to smooth out noise
       smooth_img(:,:,1) = medfilt2(img(:,:,1),[tracking_params.filter_size tracking_params.filter_size]);
       smooth_img(:,:,2) = medfilt2(img(:,:,2),[tracking_params.filter_size tracking_params.filter_size]);
       smooth_img(:,:,3) = medfilt2(img(:,:,3),[tracking_params.filter_size tracking_params.filter_size]);
       
       % for the 2nd image only
       if (i == data_params.frame_ids(1)+1)
           % check for bounds before generating search window
           y_min_tmp = (tracking_params.rect(2)+round(tracking_params.rect(4)/2)) - tracking_params.search_half_window_size;
           if (y_min_tmp < 1)
               y_min_tmp = 1;
           end
           
           y_max_tmp = (tracking_params.rect(2)+round(tracking_params.rect(4)/2)) + tracking_params.search_half_window_size;
           if y_max_tmp > size(smooth_img,1)
               y_max_tmp = size(smooth_img,1); 
           end
           
           x_min_tmp = (tracking_params.rect(1)+round(tracking_params.rect(3)/2)) - tracking_params.search_half_window_size;
           if (x_min_tmp < 1)
               x_min_tmp = 1;
           end
           
           x_max_tmp = (tracking_params.rect(1)+round(tracking_params.rect(3)/2)) + tracking_params.search_half_window_size;
           if x_max_tmp > size(smooth_img,2)
               x_max_tmp = size(smooth_img,2); 
           end
           % generating search window
           window = smooth_img(y_min_tmp:y_max_tmp, x_min_tmp:x_max_tmp,:);
       else
           % for all other images
           % check for bounds before generating search window
           y_min_tmp = (rect(2)+round(tracking_params.rect(4)/2)) - tracking_params.search_half_window_size;
           if (y_min_tmp < 1)
               y_min_tmp = 1;
           end
           
           y_max_tmp = (rect(2)+round(tracking_params.rect(4)/2)) + tracking_params.search_half_window_size;
           if y_max_tmp > size(smooth_img,1)
               y_max_tmp = size(smooth_img,1); 
           end
           
           x_min_tmp = (rect(1)+round(tracking_params.rect(3)/2)) - tracking_params.search_half_window_size;
           if (x_min_tmp < 1)
               x_min_tmp = 1;
           end
           
           x_max_tmp = (rect(1)+round(tracking_params.rect(3)/2)) + tracking_params.search_half_window_size;
           if x_max_tmp > size(smooth_img,2)
              x_max_tmp = size(smooth_img,2); 
           end
           % generating search window
           window = smooth_img(y_min_tmp:y_max_tmp, x_min_tmp:x_max_tmp,:);
       end
       figure, imshow(window)
       minxy = [];
       correlation_output = [];
       img_hist_list = [];
       for j = 1 : size(window,2)-tracking_params.rect(3)
            for k = 1 : size(window,1)-tracking_params.rect(4)
                ymax_template = k+tracking_params.rect(4)-1;
                ymin_template = k;
                % checking bound for template
                if (ymax_template > size(window,1))
                    ymax_template = size(window,1);
                    ymin_template = ymax_template - tracking_params.rect(4);
                end
                if (ymin_template < 0)
                    ymin_template = 0;
                    ymax_template = ymin_template + tracking_params.rect(4);
                end
                
                xmax_template = j+tracking_params.rect(3)-1;
                xmin_template = j;
                if (xmax_template > size(window,2))
                    xmax_template = size(window,2);
                    xmin_template = xmax_template - tracking_params.rect(3);
                end
                if (xmin_template < 0)
                    xmin_template = 0;
                    xmax_template = xmin_template + tracking_params.rect(3);
                end
                % generating template
                template = window(ymin_template:ymax_template,...
                            xmin_template:xmax_template,:);
                % weighing template with Epanechnikov Kernel
                if(size(template,1) == size(kernel,1) && size(template,2) == size(kernel,2))
                    template(:,:,1) = template(:,:,1).*kernel;
                    template(:,:,2) = template(:,:,2).*kernel;
                    template(:,:,3) = template(:,:,3).*kernel;
                end
                % generating histogram
                [ind,~] = rgb2ind(template,tracking_params.bin_n);
                img_hist = imhist(ind(:), map1);
                % calculating normalized cross correlation
                tmp_correlation_output = ncc_hist(img1_hist,img_hist);
                
                minxy = [minxy; j k];
                correlation_output = [correlation_output abs(tmp_correlation_output)];
                img_hist_list = [img_hist_list img_hist];
            end
       end
        % find the maximum correlation value and its corresponding index
        [~,idx] = max(correlation_output);
        if i == data_params.frame_ids(1)+1
            rect = [((tracking_params.rect(1)+round(tracking_params.rect(3)/2)) - tracking_params.search_half_window_size) + minxy(idx,1)...
                    ((tracking_params.rect(2)+round(tracking_params.rect(4)/2)) - tracking_params.search_half_window_size) + minxy(idx,2)...
                    tracking_params.rect(3)...
                    tracking_params.rect(4)];
        else
            rect = [((rect(1)+round(rect(3)/2)) - tracking_params.search_half_window_size) + minxy(idx,1)...
                    ((rect(2)+round(rect(4)/2)) - tracking_params.search_half_window_size) + minxy(idx,2)...
                    tracking_params.rect(3)...
                    tracking_params.rect(4)];
        end
        img = drawBox(img, rect, [0, 0, 255], 3);
        imwrite(img, fullfile(data_params.out_dir,...
        data_params.genFname(data_params.frame_ids(i))));
    end
end

%% Calculating Epanechnikov Kernel
function k = computeEpanechnikovKernel(wt,ht)
    k = zeros(ht,wt);
    for i = 1:ht
        for j = 1:wt
            tmp = ((i - ht/2)/ht)^2 + ((j - wt/2)/wt)^2;
            if tmp < 1
                k(i,j) = 1 - norm(tmp)^2;
            else
                k(i,j) = 0;
            end
        end
    end
end

%% Calculating Normalized Cross Correlation (NCC)
function d = ncc_hist(h1, h2)
    mean_h1 = mean(h1);
    mean_h2 = mean(h2);

    num = sum((h1-mean_h1).*(h2-mean_h2));
    denum = sqrt(sum((h1-mean_h1).^2))*sqrt(sum((h2-mean_h2).^2));

    d = num/denum;
end