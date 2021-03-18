function [miceMap] = miceReadMap(f)
    r=readtable(f,'ReadVariableNames',false);
    miceMap=table('Size',[size(r,1),3],'VariableTypes',{'string','cell','cell'});
    miceMap.Properties.VariableNames={'Name','N','L'};
    for i=1:size(miceMap,1)
        miceMap.Name(i)=r.Var1(i);
        if isnan(r.Var4(i)) % only one mouse
            miceMap.N{i}=r.Var2(i);
            a=strsplit(r.Var3{i},',');
            aa=zeros(length(a),1);
            for j=1:length(a)
                aa(j)=str2num(a{j});
            end
            miceMap.L{i}=aa;
        else % two mice on slide
            miceMap.N{i}=[r.Var2(i),r.Var4(i)];
            a=strsplit(r.Var3{i},',');
            b=strsplit(r.Var5{i},',');
            c=max(length(a),length(b));
            aa=zeros(c,1);
            bb=zeros(c,1);
            for j=1:length(a)
                aa(j)=str2num(a{j});
            end
            for j=1:length(b)
                bb(j)=str2num(b{j});
            end
            miceMap.L{i}=[aa,bb];
        end
    end
end

