num=[100000];
rd='X:\python\AttentionDeepMIL\image2gene\revision\run_handpicked\run2\';    % top 5

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
fprintf('rt\t%0.5f\t%0.5f\n',mean(er1),median(er1));

rd='X:\python\AttentionDeepMIL\image2gene\revision\run_handpicked_test\run2\';    % top 5
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

rd='X:\R\tb\genes\rerun\run_handpicked_random5000\';
d=dir(strcat(rd,'*_test.txt'));
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
fprintf('5\t%0.5f\t%0.5f\n',mean(er1),median(er1));



rd='X:\R\tb\genes\rerun\run_handpicked_random5000_test\';
d=dir(strcat(rd,'*_test.txt'));
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
fprintf('5\t%0.5f\t%0.5f\n',mean(er1),median(er1));



rd='X:\R\tb\genes\rerun\run_handpicked_random\';
d=dir(strcat(rd,'*_test.txt'));
n1=[];
gt1=[];
pr1=[];
for i=1:length(d)
    t=fileread(fullfile(d(i).folder,d(i).name));
    t=regexprep(t,'\n ',' ');
    t=strsplit(t,'\n');
    for j=1:length(t)-1
        tt=strsplit(t{j},'\t');
        
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
fprintf('75\t%0.5f\t%0.5f\n',mean(er1),median(er1));

rd='X:\R\tb\genes\rerun\run_handpicked_random_test\';
d=dir(strcat(rd,'*_test.txt'));
n1=[];
gt1=[];
pr1=[];
for i=1:length(d)
    t=fileread(fullfile(d(i).folder,d(i).name));
    t=regexprep(t,'\n ',' ');
    t=strsplit(t,'\n');
    for j=1:length(t)-1
        tt=strsplit(t{j},'\t');
        
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
fprintf('75\t%0.5f\t%0.5f\n',mean(er1),median(er1));