img1 = imread('ppw - 01.tif');
img1 = imresize(img1, 0.005);
img2 = imread('ppw - 02.tif');
img2 = imresize(img2, 0.005);
img3 = imread('ppw - 03.tif');
img3 = imresize(img3, 0.005);
img4 = imread('ppw - 04.tif');
img4 = imresize(img4, 0.005);
img5 = imread('ppw - 05.tif');
img5 = imresize(img5, 0.005);
img6 = imread('ppw - 06.tif');
img6 = imresize(img6, 0.005);
img7 = imread('ppw - 07.tif');
img7 = imresize(img7, 0.005);
img8 = imread('ppw - 08.tif');
img8 = imresize(img8, 0.005);
img9 = imread('ppw - 09.tif');
img9 = imresize(img9, 0.005);
img10 = imread('ppw - 10.tif');
img10 = imresize(img10, 0.005);
img11 = imread('ppw - 11.tif');
img11 = imresize(img11, 0.005);
img12 = imread('ppw - 12.tif');
img12 = imresize(img12, 0.005);
[x, y, c] = size(img1);

archive = {img1,img2,img3,img4,img5,img6,img7,img8,img9,img10,img11,img12};
iarchive = zeros(size(img1,1)*size(img1,2)*size(img1,3),12);% weighted function for img
for p =1:12
    img = archive{p};
    iarchive(:,p) = img(:);
end
for i = 1:x*y*c
B(i,:) = [log(30), log(15),log(8),log(4),log(2),log(1),log(1/2),log(1/4),log(1/8),log(1/15),log(1/30),log(1/60)];
end

zmin = 0;
zmax = 255;

%% generating weighting function
W = zeros(size(img1,1)*size(img1,2)*size(img1,3),12);% weighted function for img
w = zeros(size(img1));
for p = 1:12
    img = archive{p};
    for k = 1:c
        for i = 1:x
            for j = 1:y
                if img(i,j,k) <= 1/2*(zmin+zmax)
                    w(i,j,k) = img(i,j,k) - zmin;
                elseif img(i,j,k) > 1/2*(zmin+zmax)
                    w(i,j,k) = zmax - img(i,j,k);
                end
            end
        end
    end
    W(:,p) = w(:);
end

%% generating response funtion
G = zeros(size(img1,1)*size(img1,2)*size(img1,3),12);% weighted function for response
l=5;
G = gsolve(iarchive,B,l,W);

%% getting HDR image
R = zeros(size(img1,1)*size(img1,2)*size(img1,3),12);
D = zeros(size(img1,1)*size(img1,2)*size(img1,3),12);
for i = 1:x*y*c-1
    for j =1:12
        R(i,1) = W(i,j)*(G(i,1)-B(i,j))+ R(i,1);
        D(i,1)=  D(i,1) +  W(i,j);
    end
    R(i) = R(i)/D(i);
end


             
            
        