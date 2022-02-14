clear;clc
close all;

addpath(genpath('/set_path'))
cd('/home/data_path')

kk = 1;
flist = dir('./');
for ii= 3:length(flist)
    
    % Load the image and ROI
    fpath = strcat('/home/data_path/',flist(ii).name,'/');
    cd(fpath)    
    filelist = dir([fpath]);
    num = numel(filelist);
    temp_pat_num = flist(ii,1).name;
    pat_num = strcat(temp_pat_num);
   
    fname_img = strcat(pat_num,'_1.img');
    image = load_nii(fname_img);
    fname_roi = strcat('roi_1.img'); 
    roi = load_nii(fname_roi);

%% Image setting
    %Set the ROI [0,1]
    if roi.hdr.dime.glmin == (-1024)
        idxm = find(roi.img > (-1024));
        idxn = find(roi.img == (-1024));
        roi.img(idxm) = 1;
        roi.img(idxn) = 0;
        idx = find(roi.img > 0);
        [idx_i, idx_j, idx_k] = ind2sub(size(roi.img), idx);
        inten_list = double(image.img(idx));
    elseif roi.hdr.dime.glmin == (-2047)
        idxm = find(roi.img > (-2047));
        idxn = find(roi.img == (-2047));
        roi.img(idxm) = 1;
        roi.img(idxn) = 0;
        idx = find(roi.img > 0);
        [idx_i, idx_j, idx_k] = ind2sub(size(roi.img), idx);
        inten_list = double(image.img(idx));
    elseif roi.hdr.dime.glmin == (-4047)
        idxm = find(roi.img > (-4047));
        idxn = find(roi.img == (-4047));
        roi.img(idxm) = 1;
        roi.img(idxn) = 0;
        idx = find(roi.img > 0);
        [idx_i, idx_j, idx_k] = ind2sub(size(roi.img), idx);
        inten_list = double(image.img(idx));
    elseif roi.hdr.dime.glmin == (-1000)
        idxm = find(roi.img > (-1000));
        idxn = find(roi.img == (-1000));
        roi.img(idxm) = 1;
        roi.img(idxn) = 0;
        idx = find(roi.img > 0);
        [idx_i, idx_j, idx_k] = ind2sub(size(roi.img), idx);
        inten_list = double(image.img(idx));
    else
        idx = find(roi.img > 0);
        idxm = find(roi.img > 0);
        idxn = find(roi.img == 0);
        roi.img(idxm) = 1;
        roi.img(idxn) = 0;
        [idx_i, idx_j, idx_k] = ind2sub(size(roi.img), idx);
        inten_list = double(image.img(idx));
    end
    
    % define two sub ROIs
    se = zeros(3,3,3);  
    se(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
    se(:,:,2) = [0 1 0; 1 1 1; 0 1 0]; 
    se(:,:,3) = [0 0 0; 0 1 0; 0 0 0];
    roi_in = roi.img;
    roi_in2 = roi_in > 0;
    roi_ratio = 1;
    while roi_ratio > 2/3 % outercore has 1/3 inner core has 2/3
        roi_update = imerode(roi_in, se);
        roi_update2 = roi_update > 0;
        roi_ratio = sum(roi_update2(:))/sum(roi_in2(:));
        roi_in = roi_update;
    end
    roi_inner_core = roi_update;
    roi_outer_core = roi.img - roi_inner_core;
    
    
    % if can't calculate inner space by using 3d se (due to <3 z-slices )
    if sum(sum(sum(roi_update)))==0
        se = zeros(3,3,3);  
        se(:,:,1) = [0 0 0; 0 0 0; 0 0 0];
        se(:,:,2) = [0 0 0; 1 1 1; 0 0 0]; 
        se(:,:,3) = [0 0 0; 0 0 0; 0 0 0];
        roi_in = roi.img;
        roi_in2 = roi_in > 0;
        roi_ratio = 1;
        while roi_ratio > 2/3 % outercore has 1/3 inner core has 2/3
            roi_update = imerode(roi_in, se);
            roi_update2 = roi_update > 0;
            roi_ratio = sum(roi_update2(:))/sum(roi_in2(:));
            roi_in = roi_update;
        end
        roi_inner_core = roi_update;
        roi_outer_core = roi.img - roi_inner_core;
    end
    
    % locate ROI and compute feature, 3d
    idx = find(roi.img > 0);
    [idx_i, idx_j, idx_k] = ind2sub(size(roi.img), idx);
    inten_list = double(t1.img(idx));
    [~, ~, zz] = size(roi.img);
    temp = 0;
    for i_zz = 1:zz
        if numel(find(roi.img(:,:,i_zz) > 0)) > 0
            temp = vertcat(temp, i_zz);
        end
    end
    temp(1,:) = [];
    min_slice_zz = temp(1);
    max_slice_zz = temp(end);
    med = floor(median(min_slice_zz:max_slice_zz));
    

%% Check and save the overlay of Image and tumor ROI
    new_fpath = strcat('/home/data_path/',flist(ii).name,'/');
    cd(new_fpath)
    figure_path_pat =  new_fpath;
    
    I = img_slice;
    R = roi_slice;
    green = cat(3, zeros(size(R)),ones(size(R)), zeros(size(R))); 
    
    figure(1)
    imshow(I,'DisplayRange',[min(I(:)) max(I(:))]);
    hold on 
    h = imshow(green);
    hold off 
    set(h, 'AlphaData', R) 
    saveas(gcf, strcat('tumor_overlay.tif'));
    
    R=roi_inner_core_slice;
    figure(2)
    imshow(I,'DisplayRange',[min(I(:)) max(I(:))]);
    hold on 
    h = imshow(green);
    hold off 
    set(h, 'AlphaData', R) 
    saveas(gcf, strcat('tumor_inner.tif'));
    
    R=roi_outer_core_slice;
    figure(3)
    imshow(I,'DisplayRange',[min(I(:)) max(I(:))]);
    hold on 
    h = imshow(green);
    hold off 
    set(h, 'AlphaData', R) 
    saveas(gcf, strcat('tumor_outer.tif'));
  
    close all;
    
    %% SHAPE features
    feature_SHAPE = rad_SHAPE(roi,pat_num,figure_path_pat,'SHAPE');
    
    %% HISTOGRAM features
    % Whole HISTOGRAM features
    feature_HIST = rad_HIST(image,roi,pat_num,figure_path_pat,'CT',4096);
    
    % HISTOGRAM positive pixel features
    feature_HIST_pp = rad_HIST_pp(t1,roi,pat_num,figure_path_pat,'CT',3072);
    
    % HISTOGRAM inner/outer/delta features
    feature_HIST_diff = rad_HIST_iod(image,roi_inner_core,roi_outer_core,pat_num,figure_path_pat,'CT',4096);
    
    %% GLCM features
    % Whole GLCM
    feature_GLCM = rad_GLCM(image,roi,pat_num,figure_path_pat,'CT',128);
    
    % Whole GLCM in/out/delta features
    img2 = image.img;
    feature_GLCM_diff = rad_GLCM_iod(img2,roi_inner_core,roi_outer_core,pat_num,figure_path_pat,'CT',128);
    
    %Subsampling GLCM (Added for CT lung radiomics, 20171226)
    t1_sub = image.img(1:3:end, 1:3:end, 1:3:end);
    roi_sub = roi.img(1:3:end, 1:3:end, 1:3:end);
    feature_GLCM_sub = rad_GLCM_sub(t1_sub,roi_sub,pat_num,figure_path_pat,'CT',128);
    
    %Subsampling GLCM in/out/delta features
    roi_in_sub = roi_inner_core(1:3:end, 1:3:end, 1:3:end);
    roi_out_sub = roi_outer_core(1:3:end, 1:3:end, 1:3:end);
    feature_GLCM_sub_diff = rad_GLCM_iod_sub(t1_sub,roi_in_sub,roi_out_sub,pat_num,figure_path_pat,'CT',128);
    
    %% 3-D LoG features
    cd(fpath);
    roi_log = load_nii(fname_roi);
    [x y z] = size(roi_log.img);
    if roi_log.hdr.dime.glmin == (-1024)
        idxm = find(roi_log.img > (-1024));
        idxn = find(roi_log.img == (-1024));
        roi_log.img(idxm) = 1;
        roi_log.img(idxn) = 0;
    elseif roi_log.hdr.dime.glmin == (-2047)
        idxm = find(roi_log.img > (-2047));
        idxn = find(roi_log.img == (-2047));
        roi_log.img(idxm) = 1;
        roi_log.img(idxn) = 0;
    elseif roi_log.hdr.dime.glmin == (-4095)
        idxm = find(roi_log.img > (-4095));
        idxn = find(roi_log.img == (-4095));
        roi_log.img(idxm) = 1;
        roi_log.img(idxn) = 0;
    elseif roi_log.hdr.dime.glmin == (-1000)
        idxm = find(roi_log.img > (-1000));
        idxn = find(roi_log.img == (-1000));
        roi_log.img(idxm) = 100;
        roi_log.img(idxn) = 0;
        idx = find(roi_log.img > 0);
        [idx_i, idx_j, idx_k] = ind2sub(size(roi_log.img), idx);
        inten_list = double(image.img(idx));
    end
    
    
    temp = 0;
    for i_x = 1:x
        if numel(find(roi_log.img(i_x,:,:) > 0)) > 0
            temp = vertcat(temp, i_x);
        end
    end
    temp(1,:) = [];
    min_slice_x = temp(1);
    max_slice_x = temp(end);
    
    temp = 0;
    for i_y = 1:y
        if numel(find(roi_log.img(:,i_y,:) > 0)) > 0
            temp = vertcat(temp, i_y);
        end
    end
    temp(1,:) = [];
    min_slice_y = temp(1);
    max_slice_y = temp(end);
    
    temp = 0;
    for i_z = 1:z
        if numel(find(roi_log.img(:,:,i_z) > 0)) > 0
            temp = vertcat(temp, i_z);
        end
    end
    temp(1,:) = [];
    min_slice_z = temp(1);
    max_slice_z = temp(end);
    
    cd(fpath)
    image = load_nii(fname_img);    
    image = double(image.img) .* double(roi_log.img);
    
    log_img = image(min_slice_x:max_slice_x, min_slice_y:max_slice_y, min_slice_z:max_slice_z);
    log_roi =double(roi_log.img(min_slice_x:max_slice_x, min_slice_y:max_slice_y, min_slice_z:max_slice_z));
    
    [xx yy zz] = size(log_img);
    
    temp = 0;
    for i_zz = 1:zz
        if numel(find(log_roi(:,:,i_zz) > 0)) > 0
            temp = vertcat(temp, i_zz);
        end
    end
    temp(1,:) = [];
    min_slice_zz = temp(1);
    max_slice_zz = temp(end);
    med1 = floor(median(min_slice_zz:max_slice_zz));
   
    
    sigma = [0.5 1 1.5 2 2.5 3 3.5];
    
    Log_mean = 0;
    Log_max = 0;
    Log_median = 0;
    Log_min = 0;
    Log_entropy = 0;
    Log_uniformity = 0;
    Log_std = 0;
    Log_skewness = 0;
    Log_kurtosis = 0;
    
    for s = 1:length(sigma)
        
        siz = [5 5 5];
        sig = [sigma(s) sigma(s) sigma(s)];
        siz   = (siz-1)/2;
        
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
        
        h = h/sum(h(:));
        
        arg = (x.*x/sig(1)^4 + y.*y/sig(2)^4 + z.*z/sig(3)^4 - ...
            (1/sig(1)^2 + 1/sig(2)^2 + 1/sig(3)^2));
        h = arg.*h;
        
        h = h-sum(h(:))/prod(2*siz+1);
        
        log_img1 = imfilter(log_img,h);
        
        log_img1 = log_img1 .* log_roi;
        
        cd(new_fpath)
        figure, imshow(log_img1(:,:,med1),'InitialMagnification','fit'); colorbar;
        
        idx_log = find(log_roi > 0);
        [idx_log_i, idx_log_j, idx_log_k] = ind2sub(size(log_roi), idx_log);
        inten_list_log = double(log_img1(idx_log));
        
        Log_mean = vertcat(Log_mean, mean(inten_list_log));
        Log_max = vertcat(Log_max, max(inten_list_log));
        Log_median = vertcat(Log_median, median(inten_list_log));
        Log_min = vertcat(Log_min, min(inten_list_log));
        
        
        hist_bin_center_log = -1024:3071;
        hist_log = hist(inten_list_log(:), hist_bin_center_log);
        hist_log = hist_log./sum(hist_log);
        figure, plot( hist_bin_center_log, hist_log);
        idx_log = find(hist_log > 0);
        Log_entropy = vertcat(Log_entropy, -sum(hist_log(idx_log).*log2(hist_log(idx_log))));
        
        Log_uniformity = vertcat(Log_uniformity, sum(hist_log.^2));
        Log_std = vertcat(Log_std, std(inten_list_log));
        Log_skewness = vertcat(Log_skewness, skewness(inten_list_log));
        Log_kurtosis = vertcat(Log_kurtosis, kurtosis(inten_list_log));
    end
    close all;
    
    Log_mean(1,:) = []
    Log_max(1,:) = []
    Log_median(1,:) = []
    Log_min(1,:) = []
    Log_entropy(1,:) = []
    Log_uniformity(1,:) = []
    Log_std(1,:) = []
    Log_skewness(1,:) = []
    Log_kurtosis(1,:) = []
    
    %% 2-D features
    med1 = floor(median(min_slice_zz:max_slice_zz));
    bw = bwperim(log_img(:,:,med1));
    stats = regionprops(bw, 'Area','Perimeter', 'Eccentricity','Solidity' );
    perimeter_val = cat(1, stats.Perimeter);
    area_val = cat(1, stats.Area);
    
    %Roundness Factor
    RF = ((4*(pi))*area_val)./(perimeter_val.^2);
    
    %Eccentricity 
    Eccentricity_val = cat(1, stats.Eccentricity);
    
    %Solidity : the ratio of the tumor area over the area of the convex hull bounding the tumor.
    Solidity_val = cat(1, stats.Solidity);
    
    close all;
    
    %% 3-D fractal dimension
    cd(fpath)
    image1 = load_nii(fname_img);
    cd(fpath);
    roi_fd = load_nii(fname_roi);
    if roi_fd.hdr.dime.glmin == (-1024)
        idxm = find(roi_fd.img > (-1024));
        idxn = find(roi_fd.img == (-1024));
        roi_fd.img(idxm) = 1;
        roi_fd.img(idxn) = 0;
        
    elseif roi_fd.hdr.dime.glmin == (-2047)
        idxm = find(roi_fd.img > (-2047));
        idxn = find(roi_fd.img == (-2047));
        roi_fd.img(idxm) = 1;
        roi_fd.img(idxn) = 0;
        
    elseif roi_fd.hdr.dime.glmin == (-4095)
        idxm = find(roi_fd.img > (-4095));
        idxn = find(roi_fd.img == (-4095));
        roi_fd.img(idxm) = 1;
        roi_fd.img(idxn) = 0;
        
    elseif roi_fd.hdr.dime.glmin == (-1000)
        idxm = find(roi_fd.img > (-1000));
        idxn = find(roi_fd.img == (-1000));
        roi_fd.img(idxm) = 100;
        roi_fd.img(idxn) = 0;
        idx = find(roi_fd.img > 0);
        [idx_i, idx_j, idx_k] = ind2sub(size(roi_fd.img), idx);
        inten_list = double(image.img(idx));
    end
    image_roi = single(image1.img) .* single(roi_fd.img);
    fd_roi=image_roi(min_slice_x:max_slice_x, min_slice_y:max_slice_y, min_slice_z:max_slice_z);
    roi_fd_resize = roi_fd.img(min_slice_x:max_slice_x, min_slice_y:max_slice_y, min_slice_z:max_slice_z);
    [xx yy zz] = size(fd_roi);
    
    temp = 0;
    for i_zz = 1:zz
        if numel(find(roi_fd_resize(:,:,i_zz) > 0)) > 0
            temp = vertcat(temp, i_zz);
        end
    end
    temp(1,:) = [];
    min_slice_zz = temp(1);
    max_slice_zz = temp(end);
    
    med1 = floor(median(min_slice_zz:max_slice_zz));
    cd(new_fpath)
    figure, 
    imagesc(fd_roi(:,:,med1));colormap gray;
    saveas(gcf, strcat('roi_img_fd_',pat_num,'.tif'));
    [n, r, FD, FA, LC] = fractalanalysis(fd_roi, 'slope');
    saveas(gcf, strcat('fractal_plot',pat_num,'.tif'));
    figure, 
    imagesc(roi_fd_resize(:,:,med1));colormap gray;
    close all;
    
    
    %% Save all features to an excel file
    val = [feature_HIST,  feature_HIST_diff, feature_SHAPE, feature_GLCM, feature_GLCM_diff, feature_GLCM_sub, feature_GLCM_sub_diff, Log_mean', Log_max', Log_median', Log_min', Log_entropy', Log_uniformity', Log_std', Log_skewness', Log_kurtosis', RF(1),  Eccentricity_val(1), Solidity_val(1), FD];
    xlswrite(strcat('pat_',pat_num,'_radiomics_whole.xlsx'), val, 1, 'A2');
    kk = kk + 1;
    
end

