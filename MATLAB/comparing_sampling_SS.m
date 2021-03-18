rd='X:\R\tb\genes\revision\run2\on_predictions\';
nums=[100000];
centers=[2];

gt=csvread(strcat(rd,'gt.csv'))';
for num=nums
    fprintf('%i\n',num);
    fprintf('\tSS\t\tnSS\n');
    
    % Method 0
    t=csvread(strcat(rd,'random_',num2str(num),'.csv'));
    a=confusionmat(gt(:,2),t);
    fprintf('r\t%0.4f\t%0.4f\n',a(2,2)/sum(a(2,:)),a(1,1)/sum(a(1,:)));
    
    % Method 1
    t=csvread(strcat(rd,'uniform_',num2str(num),'.csv'));
    a=confusionmat(gt(:,2),t);
    fprintf('u\t%0.4f\t%0.4f\n',a(2,2)/sum(a(2,:)),a(1,1)/sum(a(1,:)));
    
    % Method 2
    for center=centers
        t=csvread(strcat(rd,'segments_',num2str(num),'_',num2str(center),'.csv'));
        a=confusionmat(gt(:,2),t);
        fprintf('s2\t%0.4f\t%0.4f\n',a(2,2)/sum(a(2,:)),a(1,1)/sum(a(1,:)));
    end
end


% External testing
rd='X:\R\tb\genes\revision\run2\on_predictions_test\';
nums=[100000];
centers=[2];

gt=csvread(strcat(rd,'gt.csv'))';
th=4;
for num=nums
    fprintf('%i\n',num);
    fprintf('\tSS\t\tnSS\n');
    
    % Method 0
    tgt=csvread(strcat(rd,'random_',num2str(num),'_gt.csv'));
    tpr=csvread(strcat(rd,'random_',num2str(num),'_pr.csv'));
    agt=confusionmat(gt(:,2),tgt(1,:));
    apr=confusionmat(gt(:,2),double(sum(tpr,1)>th));
    fprintf('r\t%0.4f\t%0.4f\n',apr(2,2)/sum(apr(2,:)),apr(1,1)/sum(apr(1,:)));
    fprintf('\t%0.4f\t%0.4f\n',agt(2,2)/sum(agt(2,:)),agt(1,1)/sum(agt(1,:)));
    
    % Method 1
    tgt=csvread(strcat(rd,'uniform_',num2str(num),'_gt.csv'));
    tpr=csvread(strcat(rd,'uniform_',num2str(num),'_pr.csv'));
    agt=confusionmat(gt(:,2),tgt(1,:));
    apr=confusionmat(gt(:,2),double(sum(tpr,1)>th));
    fprintf('u\t%0.4f\t%0.4f\n',apr(2,2)/sum(apr(2,:)),apr(1,1)/sum(apr(1,:)));
    fprintf('\t%0.4f\t%0.4f\n',agt(2,2)/sum(agt(2,:)),agt(1,1)/sum(agt(1,:)));
    
    % Method 2
    for center=centers
        tgt=csvread(strcat(rd,'segments_',num2str(num),'_',num2str(center),'_gt.csv'));
        tpr=csvread(strcat(rd,'segments_',num2str(num),'_',num2str(center),'_pr.csv'));
        agt=confusionmat(gt(:,2),tgt(1,:));
        agpr=confusionmat(gt(:,2),double(sum(tpr,1)>th));
        fprintf('s\t%0.4f\t%0.4f\n',apr(2,2)/sum(apr(2,:)),apr(1,1)/sum(apr(1,:)));
        fprintf('\t%0.4f\t%0.4f\n',agt(2,2)/sum(agt(2,:)),agt(1,1)/sum(agt(1,:)));
    end
end