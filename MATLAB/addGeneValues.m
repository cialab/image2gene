% Configuration
load('/isilon/datalake/cialab/scratch/cialab/tet/MATLAB/mil/ss_labels_new.mat');
slidedir='/isilon/datalake/cialab/original/cialab/image_database/d00127/Scans';
fg=0.5;
wd='/isilon/datalake/cialab/scratch/cialab/tet/MATLAB/image2gene/revision/datasets/handpicked_test/';
m=miceReadMap('../slideMap_fixed.xlsx');
%t=readtable('../data.csv');
t=readtable('../data33.csv');
gene_names=t.Properties.VariableNames(3:end);
gene_id=cellfun(@(x) str2num(x(2:end)),gene_names,'UniformOutput',false);
gene_id=cat(1,gene_id{:});
gene_id=uint32(gene_id)';
genes=[4337,8185,8187,8728,18325];

nums=[2500,10000,100000];
centers=[2];%,4,6,8];

% Loop
for n=1:size(m,1)
    mice=m.N{n};
    for i=1:length(mice)
        [genedata_idx,~]=find(mice(i)==t.MouseNum);
        idx=find(mouseNum==mice(i));
        if ~isempty(genedata_idx)
            gene_exp=table2array(t(genedata_idx,3:end));
            for num=nums
                % Method 0
                h5fn=fullfile(wd,strcat('random/',num2str(num),'/',num2str(mice(i)),'.h5'));
                h=h5info(h5fn);
                if length(h.Datasets)==1
                    h5create(h5fn,'/label',size(uint8(label)),'Datatype','uint8');
                    h5write(h5fn,'/label',uint8(label));
                    h5create(h5fn,'/gene_exp',size(gene_exp),'Datatype','double');
                    h5write(h5fn,'/gene_exp',gene_exp);
                end
                
                % Method 1
                h5fn=fullfile(wd,strcat('uniform/',num2str(num),'/',num2str(mice(i)),'.h5'));
                h=h5info(h5fn);
                if length(h.Datasets)==1
                    h5create(h5fn,'/label',size(uint8(label)),'Datatype','uint8');
                    h5write(h5fn,'/label',uint8(label));
                    h5create(h5fn,'/gene_exp',size(gene_exp),'Datatype','double');
                    h5write(h5fn,'/gene_exp',gene_exp);
                end
                
                % Method 2
                for center=centers
                    h5fn=fullfile(wd,strcat('segments/',num2str(num),'/',num2str(center),'/',num2str(mice(i)),'.h5'));
                    h=h5info(h5fn);
                    if length(h.Datasets)==1
                        h5create(h5fn,'/label',size(uint8(label)),'Datatype','uint8');
                        h5write(h5fn,'/label',uint8(label));
                        h5create(h5fn,'/gene_exp',size(gene_exp),'Datatype','double');
                        h5write(h5fn,'/gene_exp',gene_exp);
                    end
                end
            end
        end
    end
end