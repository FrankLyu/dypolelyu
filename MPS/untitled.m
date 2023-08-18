n = 2;
img = zeros(51);
for i=1:n
%     FileInfo = dir('\\BEC2-ANDOR\Users\BEC2\Desktop\dypoleImaging-fluo\andorCommand\image\combinedShots.fits');
%     if date ~= FileInfo.datenum
%         date = FileInfo.datenum;
name = append('Data-', num2str(i), '.fits');        
pic_total = fitsread(append("D:\MIT2023\Tweezer_histogram\histogram\20230724 7GHz deep tweezer on\", name));
pic_w = pic_total(:,:,1);
pic_wo = pic_total(:,:,2);

s = 150; % this is the thresold for binarization
size = 51; % size of the picture
sigma = 1; % sigma for the low pass filter
low = 0; % low and high for the colorbar
high = 0.5;
pic_w_B = pic_w>s;% binarization
pic_wo_B = pic_wo>s;
image = LowPassFilter(pic_w_B,size,sigma)-0*LowPassFilter(pic_wo_B,size,sigma); % subtract the two images and plot
image(isnan(image))=0;

% pic_final(:,:,1) = exp(image);
% % pic_final(:,:,1) = exp(image);
% pic_final(:,:,2) = ones(17);
% pic_final(:,:,3) = zeros(17);

% img=img+image;
img=img+pic_w;

pic = pic_w(21:27, 22:31)-pic_wo(21:27, 22:31) ;
sum_array(i)=sum(sum(pic));
end
%     pause(1)
% end
img=img./n;
% figure;imshow(img,[low, 0.1])
figure;imshow(img,[90, 230])
