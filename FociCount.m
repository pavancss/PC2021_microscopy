% This code performs an unbiased and unsupervised masking of DNA and CENP-A foci in an image using Otsu's algorithm.
% The DNA masks were then identified as independent DNA masses and the number of CENP-A foci in each were counted. 

files=dir(fullfile(pwd,'*.tif')); %Works only in the present working directory. 

Aggstats=struct('ImageName',{},'ChrArea',[],'ChrInt',[],'FociNumber',[],'CenInt',[],'CenArea',[]);
xctr=1;
for ctr=1:length(files)
    imstr=files(ctr).name;
    if contains(imstr,'h1dep') % String identifying each condition. Can also be simplified to count in all the images. 
        
        F=double(imread(imstr,2)); % CENP-A channel
        D=double(imread(imstr,1)); % DNA channel

        Dblur=imdilate(D,strel('disk',5));

        Dbw=imbinarize(Dblur/1e4,graythresh(Dblur/1e4));
        Dbwfilt=bwareafilter(bwlabel(Dbw),1000,100000);
        Dbg=imgavg(D.*(Dbwfilt==0));
        DbgCorr=D-Dbg*ones(size(D));
        chrmask=+(Dbwfilt>0);
        
        fbg=imgavg(F.*+(chrmask==0));
        %FbgCorr=F-Fblur;
        FbgCorr=F-fbg*ones(size(F));
        
        
        Fbw=imbinarize(FbgCorr/1e4,graythresh(FbgCorr/1e4));
        Fbwmask=Fbw.*chrmask;
        Fbwfilt=bwareafilter(bwlabel(Fbwmask),10,1000);

        imwrite(Dbwfilt,strcat(imstr(1:end-8),'DNAmask','.jpeg'),'jpeg');
        imwrite(Fbwfilt,strcat(imstr(1:end-8),'CENmask','.jpeg'),'jpeg');

        
        for chrindex=1:max(max(Dbwfilt))
            chrmask=Dbwfilt==chrindex;
            cenChr=bwlabel(Fbwfilt.*chrmask);
            Aggstats(xctr).ImageName=imstr(1:end-8);
            Aggstats(xctr).ChrArea=sum(sum(chrmask));
            Aggstats(xctr).ChrInt=imgavg(DbgCorr.*chrmask);
            Aggstats(xctr).FociNumber=max(max(cenChr));
            Aggstats(xctr).CenInt=sum(sum(FbgCorr.*(cenChr>0)));
            Aggstats(xctr).CenArea=sum(sum(cenChr>0));
            xctr=xctr+1;
        end
       
    end
    
    
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



