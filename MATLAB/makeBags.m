function makeBags(bb,b)
rng(1);
addpath('/isilon/datalake/cialab/scratch/cialab/tet/.usr/include/openslide');
addpath('/isilon/datalake/cialab/scratch/cialab/tet/.usr/lib');
openslide_load_library();

% Configuration
load('/isilon/datalake/cialab/scratch/cialab/tet/MATLAB/mil/ss_labels_new.mat');
slidedir='/isilon/datalake/cialab/original/cialab/image_database/d00127/Scans';
wd='/isilon/datalake/cialab/scratch/cialab/tet/MATLAB/image2gene/revision/datasets/224_5percent'; %bySlide
m=miceReadMap('../slideMap_fixed.xlsx');
avoid=[194,196,198,200,202,204,206,208,211,214,216,218,221,224,231,240,250,252,256,261,272,288,299,305,311,350,355,368,375,384,394,407,411,418,421,426,440,455,474,195,197,199,201,203,205,207,210,213,215,217,220,222,227,239,249,251,253,260,265,285,296,300,307,344,351,359,372,383,386,398,409,417,419,424,429,443,471,209,225,229,232,233,234,236,241,267,274,279,281,286,295,303,316,318,329,336,340,365,367,378,389,415,422,434,437,451,452,459,465,467];
prop=0.05;

% Loop
m=m(randperm(size(m,1)),:);
l=round(linspace(0,size(m,1),bb+1));
for n=l(b)+1:l(b+1)
    slide=m.Name(n);
    mice=m.N{n};
    
    % Open file
    fp=openslide_open(fullfile(slidedir,char(slide)));
    [y,x]=openslide_get_level0_dimensions(fp);
    
    for i=1:length(mice)
        if ~isempty(find(mouseNum==mice(i)))
            if sum(mice(i)==avoid)==0
                fprintf('%s: %i\n',char(slide),mice(i));
                tic;

                % Da mask
                mask=imread(fullfile('../mouseMasks',strcat(num2str(mice(i)),'.png')));

                % Da bag
                makeBag_grid_sampling(fp,mice(i),y,x,mask,prop,224,224,wd);
                
                toc;
            end
        end
    end
    openslide_close(fp);
end

end