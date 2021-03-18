nums=[10000,20000,50000,100000];
centers=[2,4,6,8];
nums=[2500,10000,100000];
centers=[2];
rd='X:\python\AttentionDeepMIL\image2gene\revision\run_handpicked\run2\';    % top 5
rd='X:\python\AttentionDeepMIL\image2gene\revision\run_handpicked2\';   % random 10 genes
rd='X:\python\AttentionDeepMIL\image2gene\revision\run_handpicked3\run1\';   % random 1000 genes
%rd='X:\python\AttentionDeepMIL\image2gene\revision\run_handpicked_test\run2\';    % top 5
%rd='X:\python\AttentionDeepMIL\image2gene\revision\run_handpicked2_test\';   % random 10 genes
rd='X:\python\AttentionDeepMIL\image2gene\revision\run_handpicked3_test\run1\';   % random 1000 genes
for num=nums
    fprintf('%i\n',num);
    fprintf('\tMean\tMedian\n');
    % Method 0
    d=dir(strcat(rd,'random_',num2str(num),'_*_test.txt'));
    n1=[];
    gt1=[];
    pr1=[];
    for i=1:length(d)
        t=fileread(fullfile(d(i).folder,d(i).name));
        t=regexprep(t,'\n ',' ');
        t=strsplit(t,'\n');
        for j=1:length(t)-1
            tt=strsplit(t{j},'\t');
            if strcmp(tt{1},'collision. trying again.')
                continue
            end
            n1=cat(1,n1,str2num(tt{1}));
            gt1=cat(1,gt1,str2num(tt{2}));
            pr1=cat(1,pr1,str2num(tt{3}));
        end
    end
    er1=zeros(size(gt1,1),1);
    for i=1:length(er1)
        er1(i)=sum(gt1(i,:).*pr1(i,:))/(norm(gt1(i,:))*norm(pr1(i,:)));
        %er1(i)=mean(abs(gt1(i,:)-pr1(i,:)).^2);
    end
    fprintf('r\t%0.5f\t%0.5f\n',mean(er1),median(er1));
    
    % Method 1
    d=dir(strcat(rd,'uniform_',num2str(num),'_*_test.txt'));
    n2=[];
    gt2=[];
    pr2=[];
    for i=1:length(d)
        t=fileread(fullfile(d(i).folder,d(i).name));
        t=regexprep(t,'\n ',' ');
        t=strsplit(t,'\n');
        for j=1:length(t)-1
            tt=strsplit(t{j},'\t');
            if strcmp(tt{1},'collision. trying again.')
                continue
            end
            n2=cat(1,n2,str2num(tt{1}));
            gt2=cat(1,gt2,str2num(tt{2}));
            pr2=cat(1,pr2,str2num(tt{3}));
        end
    end
    er2=zeros(size(gt2,1),1);
    for i=1:length(er2)
        er2(i)=sum(gt2(i,:).*pr2(i,:))/(norm(gt2(i,:))*norm(pr2(i,:)));
        %er2(i)=mean(abs(gt2(i,:)-pr2(i,:)).^2);
    end
    fprintf('u\t%0.5f\t%0.5f\n',mean(er2),median(er2));
    
    % Method 2
    for center=centers
        d=dir(strcat(rd,'segments_',num2str(num),'_',num2str(center),'_*_test.txt'));
        n3=[];
        gt3=[];
        pr3=[];
        for i=1:length(d)
            t=fileread(fullfile(d(i).folder,d(i).name));
            t=regexprep(t,'\n ',' ');
            t=strsplit(t,'\n');
            for j=1:length(t)-1
                tt=strsplit(t{j},'\t');
                if strcmp(tt{1},'collision. trying again.')
                    continue
                end
                n3=cat(1,n3,str2num(tt{1}));
                gt3=cat(1,gt3,str2num(tt{2}));
                pr3=cat(1,pr3,str2num(tt{3}));
            end
        end
        er3=zeros(size(gt3,1),1);
        for i=1:length(er3)
            er3(i)=sum(gt3(i,:).*pr3(i,:))/(norm(gt3(i,:))*norm(pr3(i,:)));
            %er3(i)=mean(abs(gt3(i,:)-pr3(i,:)).^2);
        end
        fprintf('s%i\t%0.5f\t%0.5f\n',center,mean(er3),median(er3));
        fprintf('vs_r\tp=%0.5f\n',ranksum(er3,er1,'tail','right'));
        fprintf('vs_u\tp=%0.5f\n',ranksum(er3,er2,'tail','right'));
    end
end