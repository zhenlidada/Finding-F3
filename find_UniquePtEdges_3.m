function UniquePtEdges_new = find_UniquePtEdges_3(NodeCount,UniquePtEdges,UniqueNodes)
% This is function is used to fix those points connnected with muitiple points

ExtraNodeidx = find(NodeCount==3);
inpoint_all = zeros(size(ExtraNodeidx,2),4);
size_p = zeros(size(ExtraNodeidx,2),4);

for p = 1:size(ExtraNodeidx,2)
    
    F1_4 = UniquePtEdges(find(sum(UniquePtEdges==ExtraNodeidx(p),2)),:);
    F1_4 = F1_4(F1_4~=ExtraNodeidx(p));
    
    for i =1:2
        nextone = F1_4(i);
        ran1_all = [];
        ran1_all = ExtraNodeidx(p);
        ran1_all = [ran1_all,nextone(1)];
        k=2;
        
        while ~isempty(nextone)
            ran1_all(k,1) = nextone;
            nextone = UniquePtEdges(find(sum(UniquePtEdges==nextone,2)),:);
            nextone = nextone(nextone~=ran1_all(k,1)& ~ismember(nextone,ran1_all));
            if ~isempty(nextone)
                if size(nextone,1)~=1
                    nextone = UniquePtEdges(find(sum(UniquePtEdges==ran1_all(end,1),2)),:);
                    nextone = intersect(nextone(1,:),nextone(2,:));
                    if nextone == ran1_all(1)
                        ran1_all(k,:)=[];
                        break
                    else
                        if isempty(find(ExtraNodeidx==nextone))
                            ran1_all(k,:)=[];
                            inpoint_all(p,i) = 10;
                            break
                        else
                            indt = find(ExtraNodeidx==nextone);
                            ran1_all(k,:)=[];
                            inpoint_all(p,i) = indt;
                            break
                        end
                    end
                else
                    ran1_all(k,2)=nextone;
                    k=k+1;
                    ran1_all(k,1)=nextone;
                end
            else
                ran1_all(k,2)=ran1_all(1);
                nextone = [];
            end
        end
        cell_all{p,1}{i} = ran1_all;
        size_p(p,i) = size(cell_all{p,1}{i},1);
    end
    
end

if max(size_p(:))< size(UniquePtEdges,1)/3
    indlast_end=[];
    
    for p = 1:size(ExtraNodeidx,2)
        for i =1:2
            indlast_end = [indlast_end;cell_all{p,1}{1,i}];
        end
    end
    
    indlast_end = unique(indlast_end);
    indlast_remo = [];
    for j=1:size(indlast_end,1)
        indlast_remo = [indlast_remo;find(sum(UniquePtEdges==indlast_end(j),2)>=1)];
    end
    
    UniquePtEdges(indlast_remo,:)=[];
    UniquePtEdges_new = UniquePtEdges;
    if size(unique(UniquePtEdges_new),1)~=size(UniquePtEdges_new,1)
        UniquePtEdges_new = add_edge(UniquePtEdges_new);
    end
    
else
    
    [point,edge] = find(size_p==max(size_p(:)));
    P1 = point(1);
    edge1 = edge(1);
    
    if isempty(find(inpoint_all(P1,edge1)))
        UniquePtEdges_new = cell_all{P1,1}{1,edge1};
    else
        UniquePtEdges_new = cell_all{P1,1}{1,edge1};
        P2 = inpoint_all(P1,edge1);
        P_P1 = inpoint_all(P1,:);
        P_P2 = inpoint_all(P2,:);
        P_P2(find(size_p(P2,:)==max(size_p(P2,:))))=0;
        P1_P2 = intersect(P_P1,P_P2);
        P1_P2(P1_P2==0)=[];
        P3 = P1_P2;
        
        if ~isempty(P3)
            E_P2_P3 = find(inpoint_all(P2,:)==P3);
            UniquePtEdges_new23 = cell_all{P2,1}{1,E_P2_P3};
            E_P1_P3 = find(inpoint_all(P1,:)==P3);
            UniquePtEdges_new3 = cell_all{P1,1}{1,E_P1_P3};
            UniquePtEdges_new3 = flip(flip(UniquePtEdges_new3,1),2);
        else
            E_P1_P2 = find(inpoint_all(P1,:)==P2);
            E_P1_P2 = E_P1_P2(E_P1_P2~=edge1);
            UniquePtEdges_new3 = cell_all{P1,1}{1,E_P1_P2};
            UniquePtEdges_new3 = flip(flip(UniquePtEdges_new3,1),2);
        end
        
        if ~isempty(P3)
            UniquePtEdges_new = [UniquePtEdges_new;UniquePtEdges_new23;UniquePtEdges_new3];
        else
            UniquePtEdges_new = [UniquePtEdges_new;UniquePtEdges_new3];
        end
    end
end

return

