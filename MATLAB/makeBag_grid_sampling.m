function makeBag_grid_sampling(fp,fn,y0,x0,mask,prop,ps,stride,wd)        
% Save time by computing boundaries for query points
boundaries=bwboundaries(mask);
minx=size(mask,2);miny=size(mask,1);maxx=0;maxy=0;
for i=1:length(boundaries)
    b=boundaries{i};
    if min(b(:,2))<minx
        minx=min(b(:,2));
    end
    if min(b(:,1))<miny
        miny=min(b(:,1));
    end
    if max(b(:,2))>maxx
        maxx=max(b(:,2));
    end
    if max(b(:,1))>maxy
        maxy=max(b(:,1));
    end
end
psx=ps*size(mask,1)/double(x0);
psy=ps*size(mask,2)/double(y0);
stride=stride*size(mask,1)/double(x0);

% [y2,x2]=openslide_get_level_dimensions(fp,2);
% im2=openslide_read_region(fp,0,0,y2,x2,'level',2);
% im2=im2(:,:,2:4);
% figure;
% imshow(im2);
% hold on;

% Grid of points
ys=miny:stride:maxy;
xs=minx:stride:maxx;
[ys,xs]=meshgrid(ys,xs);
pts=cat(2,ys(:),xs(:));
r=randsample(1:size(pts,1),size(pts,1)-round(prop*size(pts,1)));
pts(r,:)=[];

% Here's where we put things
bag=zeros(224,224,3,size(pts,1),'uint8');
keep=zeros(size(pts,1),1);
for p=1:size(pts,1)
    % Query points
    rx=pts(p,2);
    ry=pts(p,1);

    % Checks if inside the mask image
    if ((size(mask,2)>rx+psx) && (size(mask,1)>ry+psy)) && ...
        mask(round(ry),round(rx))==1 && ...
        mask(round(ry+psy),round(rx))==1 && ...
        mask(round(ry),round(rx+psx))==1 && ...
        mask(round(ry+psy),round(rx+psx))==1
        im=openslide_read_region(fp,round(rx*double(y0)/size(mask,2)),round(ry*double(x0)/size(mask,1)),ps,ps);
        im=im(:,:,2:4);
        bag(:,:,:,p)=im;
        keep(p)=1;
    end
end
bag(:,:,:,~keep)=[];
pts(~keep,:)=[];
pts(:,2)=round(pts(:,2).*(double(y0)/size(mask,2)));
pts(:,1)=round(pts(:,1).*(double(x0)/size(mask,1)));
pts=uint32(pts);

h5fn=fullfile(wd,strcat(num2str(fn),'.h5'));
h5create(h5fn,'/bag',size(bag),'Datatype','uint8');
h5write(h5fn,'/bag',bag);
h5create(h5fn,'/coords',size(pts),'Datatype','uint32');
h5write(h5fn,'/coords',pts);

end