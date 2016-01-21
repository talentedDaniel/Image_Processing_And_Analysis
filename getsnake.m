
function [xs, ys] = getsnake(image)
% image: The image on which snake is to be initialized
% xs, xy: Initial position of the snake

hold on; %Hold the image steady       
xs = [];
ys = [];

% xold = 0;
% yold = 0;

but = 1;
hold on
% Initially, the list of points is empty.
xy = [];
n = 0;
% Loop, picking up the points.
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
while but == 1
    [xi,yi,but] = ginput(1); %pick a new input
    plot(xi,yi,'ro')
    n = n+1;
     % xy is 2*n matrix.
    xy(:,n) = [xi;yi];
end
% close snake
n = n+1;
xy(:,n) = [xy(1,1);xy(2,1)];

% Interpolate with a spline curve and finer spacing.
t = 1:n;
ts = 1: 0.1: n;
xys = spline(t,xy,ts);

xs = xys(1,:);
ys = xys(2,:);

plot(xy(1,:),xy(2,:),'or',xys(1,:),xys(2,:),'-b'), axis equal

