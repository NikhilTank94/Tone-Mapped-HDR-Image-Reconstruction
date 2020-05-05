clear all
img1 = imread('ppw - 01.jpg');
img2 = imread('ppw - 02.jpg');
img3 = imread('ppw - 03.jpg');
img4 = imread('ppw - 04.jpg');
img5 = imread('ppw - 05.jpg');
img6 = imread('ppw - 06.jpg');
[x, y, c] = size(img1);
archive = {img1,img2,img3,img4,img5,img6};
B = [log(30); log(15);log(8);log(4);log(2);log(1)];

%% generating weighting function
zmin = 0;
zmax = 255;
w = zeros(256,1);
for i = 0:255
    if i <= 1/2*(zmin+zmax)
       w(i+1) = i - zmin;
    elseif i > 1/2*(zmin+zmax)
       w(i+1) = zmax - i;
    end
end  

%% generating response funtion
inew=zeros(256,6);
for z =1:c
    for p =1:6
        arc_temp = archive{p}(:,:,z);
        inew(:,p) = arc_temp(1:256);   
    end
     g = gsolve(inew,B,5,w);
     G(:,z) = g;
end

%% getting HDR image
N = zeros(x,y,c);
D = zeros(x,y,c);
I = zeros(x,y,c);
for z = 1:c
    for i = 1:x
        for j = 1:y
            for p =1:6
                N(i,j,z) = w(archive{p}(i,j,z)+1)*(G(archive{p}(i,j,z)+1,z) - B(p))+ N(i,j,z);
                D(i,j,z)=  w(archive{p}(i,j,z)+1) +  D(i,j,z);
            end
            I(i,j,z) = exp(N(i,j,z)/D(i,j,z));
        end
    end
end
I = tonemap(I);
imwrite(I, 'HDR_ghost.jpg');


