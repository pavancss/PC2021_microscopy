%{
Main function for going through a folder to analyze the images. Analysis is
only on the highest contrast z-slice of the stack. 
 
Image folder should only contain tif files to be analyzed. Tif files should be a
arranged with all channels of each z-section followed by the next (Fiji
output is good). 

Output- files is a directory with the file names in one of the columns. The
image intensities and some normalizations are in the IFData array. 

%}

files=dir(fullfile(pwd,'*.tif'));

IFData=zeros(length(files),7);
totC=4; % Number of channels in the image

for ctr=1:length(files)
    imstr=files(ctr).name;
    
    z=Zselect(imstr,totC,1); % 1 is the channel to perform the contrast selection on. Change as needed. 
    
    D=double(imread(imstr,(z-1)*totC+1));
    F=double(imread(imstr,(z-1)*totC+2));
    T=double(imread(imstr,(z-1)*totC+3));
    C=double(imread(imstr,(z-1)*totC+4));
    
    
    maskBG=+(D<600); % Select this number so that it's in between background and object intensities. 
    Dbg=mean(mean(D.*maskBG));
    Fbg=mean(mean(F.*maskBG));
    Tbg=mean(mean(T.*maskBG));
    Cbg=mean(mean(C.*maskBG));
    
    D=D-Dbg*ones(size(D));
    F=F-Fbg*ones(size(D));
    T=T-Tbg*ones(size(D));
    C=C-Cbg*ones(size(D));
    
    mask=+(D>1000); % I decide this after looking at images and set a single number for the entire experiment.
   
    IFData(ctr,1)=imgavg(D.*mask);
    IFData(ctr,2)=imgavg(T.*mask);
    IFData(ctr,3)=imgavg(F.*mask);
    IFData(ctr,4)=imgavg(C.*mask);
    IFData(ctr,5)=IFData(ctr,3)/IFData(ctr,2);
    IFData(ctr,6)=IFData(ctr,4)/IFData(ctr,2);
    IFData(ctr,7)=IFData(ctr,3)/IFData(ctr,1);
end


% This is a function to calculate the mean intensity of a 2-D image of any
% size. The function ignores all zero value pixels unlike the mean
% function. 
% Input- any 1 or 2-D array with positive values to get the average of
% Output- Mean of non-zero values in input.


function avg=imgavg(I)

[a,b]=size(I);

s=0;
n=0;

binI=+(I>0);
I=I.*binI;
for i=1:a
    s=s+sum(I(i,1:b));
    n=n+sum(binI(i,1:b));
end
%avg=s;
avg=s/n;

end

%{ 
Function to select highest contrast slice of a z-stack.
Inputs- imstr is the name of the image. totC is the number of channels in
the image. c is the channel for which the maximum projection needs to be
performed. 
Output- z value of highest contrast of channel 'c' of image imstr.
%}


function zsel=Zselect(imstr,totC,c)

zsize=numel(imfinfo(imstr))/totC;
S=zeros(zsize,1);
for ctr=1:zsize 
    F=double(imread(imstr,(ctr-1)*totC+c));
    F=F.*(F>700); % General mask to remove most background. Keep high to get actual object pixels for max intensity. 
    [r,c,val]=find(F);
    S(ctr)=mean(val);
end
[test,zsel]=max(S);


end


    
    