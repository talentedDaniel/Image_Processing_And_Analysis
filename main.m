 img=imread('test_image\coins.tif');
 imshow(img);
 [xs, ys] = getsnake(img);
 sigma =1.00;
 img = filter_function(img, sigma);
 alpha =0.40;
 beta =0.20;
 gamma =1.00;
 weline =0.30;
 weedge =0.40;
 weterm =0.70;
 inter =300;
 snake_algorithm(img, xs, ys, alpha, beta, gamma, weline, weedge, weterm, inter);

 
 
