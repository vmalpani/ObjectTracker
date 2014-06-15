function result = computeFlow(img1, img2, win_radius, template_radius, grid_MN)
    img1 = im2double(img1);
    img2 = im2double(img2);
    
    % median filtering to smooth out noise
    smooth_img1 = medfilt2(img1,[5 5]);
    smooth_img2 = medfilt2(img2,[5 5]);
    
    % normalizing image
    norm_img1 = (smooth_img1-mean(smooth_img1(:))) ./ sqrt(var(smooth_img1(:)));
    norm_img2 = (smooth_img2-mean(smooth_img2(:))) ./ sqrt(var(smooth_img2(:)));
    
    % padding only image 2 so that template fits in
    padded_img2 = padarray(norm_img2,[win_radius win_radius], 'replicate');
    fh1 = figure;
    imshow(img1), hold on
    
    % for sparse normalized cross correlation
    % plotting optical flow only on the specified grid size
    step = round(size(img1)./grid_MN);
    count = 0;
    
    x = zeros(1,size(norm_img1,1)*size(norm_img1,2)/(step(1,1)*step(1,2)));
    y = zeros(1,size(norm_img1,1)*size(norm_img1,2)/(step(1,1)*step(1,2)));
    u = zeros(1,size(norm_img1,1)*size(norm_img1,2)/(step(1,1)*step(1,2)));
    v = zeros(1,size(norm_img1,1)*size(norm_img1,2)/(step(1,1)*step(1,2)));
    
    for i = template_radius+1 : step(1,2) : size(norm_img1,2)-template_radius-1
       for j = template_radius+1 : step(1,1) : size(norm_img1,1)-template_radius-1
           count = count+1;
           
%            uncomment to plot bounds of template           
%            plot(i-template_radius,j-template_radius,'r*');
%            plot(i+template_radius,j+template_radius,'g*');

           template = norm_img1(j-template_radius:j+template_radius,i-template_radius :i+template_radius);
           window = padded_img2(j:j+2*win_radius,i:i+2*win_radius);
           
           correlation_output = normxcorr2(template, window);
           
           [~,idx] = max(correlation_output(:));
           [max_y, max_x] = ind2sub(size(correlation_output),idx);

           % (max_x - template_radius -> max in template frame
           % (i + max_x - template_radius) -> max in padded image frame
           % (i + max_x - template_radius) - (win_radius + 1) -> max in original image frame
           % max - i, max - j -> optical flow vectors
           u(count) = (i + max_x - template_radius) - (win_radius + 1) - i;
           v(count) = (j + max_y - template_radius) - (win_radius + 1) - j;
           
           x(count) = i;
           y(count) = j;
       end
    end
    quiver(x,y,u,v,0,'r');
    annotated_img = saveAnnotatedImg(fh1);
    result = annotated_img;
end

%%
function annotated_img = saveAnnotatedImg(fh)
figure(fh); % Shift the focus back to the figure fh

% The figure needs to be undocked
set(fh, 'WindowStyle', 'normal');

% The following two lines just to make the figure true size to the
% displayed image. The reason will become clear later.
img = getimage(fh);
truesize(fh, [size(img, 1), size(img, 2)]);

% getframe does a screen capture of the figure window, as a result, the
% displayed figure has to be in true size. 
frame = getframe(fh);
frame = getframe(fh);
pause(0.5); 
% Because getframe tries to perform a screen capture. it somehow 
% has some platform depend issues. we should calling
% getframe twice in a row and adding a pause afterwards make getframe work
% as expected. This is just a walkaround. 
annotated_img = frame.cdata;
end