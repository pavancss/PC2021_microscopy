% This code calculates the surface area and volume of a 3D object from a z-stack
% The starting image needs to be of the whole object from a confocal stack or a 
% deconvolved image.
% You can chose which channel to perform the segmentation on. 


files=dir(fullfile(pwd,'*.tif')); % This runs in the present directory only.

totC=3; %Total number of channels in image.
Filteredstats=struct('ImageName',{},'Volume',[],'IntegratedInt',[],'SurfaceArea',[]);
xctr=1;

for ctr=1:length(files)
    imstr=files(ctr).name;
    Aint=Interpolate3(imstr,totC,1); % Interpolate the z-stack to better estimate the volume
    Abg=pctl(Aint,95); % Background for subtraction
    Asub=Aint-Abg*ones(size(Aint)); 
    Abin=imbinarize(Asub/1e6,graythresh(Asub/1e6));  % Segment in 3D using Otsu's algorithm. Normalization required for better operation.
    
    Abgsub=Aint-pctl(Aint,5)*ones(size(Aint));
    Amasked=Abgsub.*+(Abgsub>0);
    
    stats=regionprops3(Abin,'Volume','SurfaceArea','VoxelList');
    [a,b]=size(stats);
    %label=bwlabeln(Abin);
    for regctr=1:a
        if stats.Volume(regctr)>5000
            
            Filteredstats(xctr).ImageName=imstr(1:end-10);
            Filteredstats(xctr).Volume=stats.Volume(regctr);
            Filteredstats(xctr).SurfaceArea=stats.SurfaceArea(regctr);
            
            voxellist=cell2mat(stats.VoxelList(regctr));
            mask=zeros(size(Aint));
            for voxctr=1:length(voxellist)
                mask(voxellist(voxctr,2),voxellist(voxctr,1),voxellist(voxctr,3))=1;
            end
            Filteredstats(xctr).IntegratedInt=sum(Amasked.*mask,'all');
            
            imsave_3D(rot90(Amasked.*+mask),strcat(imstr(1:end-12),'_',num2str(regctr),'.tif'),'C:\Users\Pavan\Google Drive\Funabiki-Data\Lab notebook\Analysis\B3-110\'); % Stores segmented image with interpolated intensities. 
            xctr=xctr+1;
        end
    end
 
end

% Writing percentile function for only upto 2D arrays
function xptl=pctl(A,x)
    [a,b]=size(A);
    if a>1 && b>1
        A=reshape(A,[a*b,1]);
    end
    Asort=sort(A);
    xptl=Asort(round(a*b*x/100),1);

end

% 3D Interpolation function. This was useful to generate a more continuous object than % generating big trapezoids. Can specify which channel number to perform on with c 
% input. Output is interpolated stack of just the single channel. 

function Aint=Interpolate3(imstr,totC,c)
info=imfinfo(imstr);
row=info(1).Height;
col=info(1).Width;
zsize=numel(info)/totC;

A=zeros(row,col,zsize);
for ctr=1:zsize
    A(:,:,ctr)=double(imread(imstr,(ctr-1)*totC+c));
end
intsize=2*(zsize-1)+zsize;
Aint=zeros(row,col,intsize);
for ctr=1:intsize
    
    if mod(ctr,3)==1
        z=(ctr-1)/3+1;
        Aint(:,:,ctr)=A(:,:,z);
    elseif mod(ctr,3)==2
        z=(ctr-2)/3+1;
        Aint(:,:,ctr)=(2*A(:,:,z)+A(:,:,z+1))/3;
    else
        Aint(:,:,ctr)=(A(:,:,ctr/3)+2*A(:,:,ctr/3+1))/3;
    end
       
end

end

% This is a script to save 3-D matrices as a multi page tif stack.
% Inputs- A is a 3D matrix where the 3rd dimension is the stacking
% dimension
% filename is the name of the file with appendage .tif
% dirstr is an optional argument if the files need to be saved in a
% different directory instead.

function imsave_3D(A,filename,dirstr)
[a,b,c]=size(A);
if nargin>2
   fname=strcat(dirstr,filename);
else
    fname=strcat(filename);
end

imwrite(uint16(A(:,:,1)),fname,'tif');
if c>1
    for ctr=2:c
        imwrite(uint16(A(:,:,ctr)),fname,'tif','WriteMode','append');
    end
end

end

% Takes bwlabel input and area cutoff and pares down the segmented image to
% exclude small things.

% !!!!!!!! A has to be bwlabel image. 



function Acut=bwareafilter(A,lowcutoff,hicutoff)
    
list=unique(A);
Acut=A;
for ctr=2:length(list)
    mask=A==list(ctr);
    size=sum(sum(mask));
    if size<lowcutoff
        Acut(mask)=0;
    elseif size>hicutoff
        Acut(mask)=0;
    end
end

Acut=bwlabel(Acut);

end


    

 