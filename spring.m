 image=imread('test_image\circle.jpg');
 figure;imshow(image);
 [xs, ys] = getsnake(image);
 sigma =1.00;
 image = filter_function(image, sigma);
 alpha =1;
 beta =1;
 gamma =1.00;
 wl =0.30;
 we =0.40;
 wt =0.70;
 interation =100;

 anchor_x =80;
 anchor_y =80;
 plot(anchor_x, anchor_y, 'or');
 
 [n,m]= size(xs);
 index =1;
 length = (anchor_x-xs(1,1))*(anchor_x-xs(1,1))+(anchor_y-ys(1,1))*(anchor_y-ys(1,1));
 
 for i=1:m
    p = [xs(1,i),ys(1,i)];
    temp = (anchor_x-xs(1,i))*(anchor_x-xs(1,i))+(anchor_y-ys(1,i))*(anchor_y-ys(1,i));
    if(temp<length)
        index=i;
        length =temp;
    end  
 end
 
% image: This is the image data
% xs, ys: The initial snake coordinates
% alpha: Controls tension
% beta: Controls rigidity
% gamma: Step size
% kappa: Controls enegry term
% wl, we, wt: Weights for line, edge and terminal enegy components
% iterations: No. of iteration for which snake is to be moved


%parameters
pic = image;
N = interation;
% Calculating size of image
[row col] = size(image);

%Computing external forces

eline = pic; %eline is the image intensities

[grady,gradx] = gradient(pic);
eedge = -1 * sqrt ((gradx .* gradx + grady .* grady)); %eedge is measured by gradient in the image

%masks for taking various derivatives
m1 = [-1 1];
m2 = [-1;1];
m3 = [1 -2 1];
m4 = [1;-2;1];
m5 = [1 -1;-1 1];

cx = conv2(pic,m1,'same');
cy = conv2(pic,m2,'same');
cxx = conv2(pic,m3,'same');
cyy = conv2(pic,m4,'same');
cxy = conv2(pic,m5,'same');

for i = 1:row
    for j= 1:col
        eterm(i,j) = (cyy(i,j)*cx(i,j)*cx(i,j) -2 *cxy(i,j)*cx(i,j)*cy(i,j) + cxx(i,j)*cy(i,j)*cy(i,j))/((1+cx(i,j)*cx(i,j) + cy(i,j)*cy(i,j))^1.5);
    end
end

% imview(eterm);
% imview(abs(eedge));
eext = (wl*eline + we*eedge - wt * eterm); %eext as a weighted sum of eline, eedge and eterm

[fx, fy] = gradient(eext); %computing the gradient

%initializing the snake
xs=xs';
ys=ys';
[m n] = size(xs);%size of snake
[mm nn] = size(fx);%size of original image
    
%populating the penta diagonal matrix
A = zeros(m,m);  %A matrix
D1=zeros(m,m);   %first derivative matrix
D2=zeros(m,m);   %second derivative matrix
A_alpha = zeros(m,m);  %A? matrix
A_beta = zeros(m,m);   %A? matrix
brow1=zeros(1,m);
brow1(1,1)=[1];         
brow1(1,m)=-1;        % populating a template row for first derivative matrix
brow2=zeros(1,m);
brow2(1,1:2)=[-2,1];
brow2(1,m)=1;         % populating a template row for second derivative matrix
brow_alpha=zeros(1,m);
brow_alpha(1,1:2)=[2,-1];
brow_alpha(1,m)=-1;    % populating a template row for A? matrix
brow_beta=zeros(1,m);
brow_beta(1,1:3)=[6,-4,1];
brow_beta(1,m-1:m)=[1,-4];     % populating a template row for A? matrix
brow_A=zeros(1,m);
brow_A=alpha*brow_alpha+beta*brow_beta;   % populating a template row for A matrix


for i=1:m
    D1(i,:)=brow1;
    brow1=circshift(brow1',1)';
    
    D2(i,:)=brow2;
    brow2=circshift(brow2',1)';
    
    A_alpha(i,:)=brow_alpha;
    brow_alpha=circshift(brow_alpha',1)';
    
    A_beta(i,:)=brow_beta;
    brow_beta=circshift(brow_beta',1)';
    
    A(i,:) = brow_A;
    brow_A = circshift(brow_A',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
end

[L U] = lu(A + gamma .* eye(m,m));
Ainv = inv(U) * inv(L); % Computing Ainv using LU factorization
 
%moving the snake
%VE=zeros(1,N);
%VEint=zeros(1,N);


for i=1:N;
    %Eint=(alpha*(norm(D1*xs)^2+norm(D1*ys)^2)+beta*(norm(D2*xs)^2+norm(D2*ys)^2))/2;
    %Eext=sum(interp2(eext,xs,ys));
    %VE(1,i)=Eint+Eext;
    %VEint(1,i)=Eint;
    ssx = gamma*xs - 0.05*interp2(fx,xs,ys)+0.05*(anchor_x-xs(index,1))*(anchor_x-xs(index,1));
    ssy = gamma*ys - 0.05*interp2(fy,xs,ys)+0.05*(anchor_y-ys(index,1))*(anchor_y-ys(index,1));
    %ssx = gamma*xs ;
    %ssy = gamma*ys ;
    
    %calculating the new position of snake
    xs = Ainv * ssx;
    ys = Ainv * ssy;
    
    if mod (i,25)==0 && i<1000
    imshow(image,[]); 
    hold on;
    plot([xs; xs(1)], [ys; ys(1)], '-b');
    plot(anchor_x, anchor_y, 'or');
    hold off;
	end
    pause(0.001)    
end;
%figure(1); plot(1:N,VE,'r',1:N,VEint,'g',1:N,VE-VEint,'b');
%set(gca,'fontsize',15)
%xlabel('Iterations'); ylabel('Energy'); 
%legend('Total Energy','Internal Energy','External Energy');

 
 
