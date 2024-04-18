function UniquePtEdges_new = find_UniquePtEdges_4(ExtraNodeidx,UniquePtEdges,UniqueNodes)
% This is function is used to fix those points connnected with muitiple points

inpoint_all = zeros(size(ExtraNodeidx,2),4);
size_p = zeros(size(ExtraNodeidx,2),4);

for p = 1:size(ExtraNodeidx,2)
    
    F1_4 = UniquePtEdges(find(sum(UniquePtEdges==ExtraNodeidx(p),2)),:);
    F1_4 = F1_4(F1_4~=ExtraNodeidx(p));
    
    for i =1:4
        
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

if max(size_p(:))< size(UniquePtEdges,1)/3 %devide 2 to 3
    indlast_end=[];
    
    for p = 1:size(ExtraNodeidx,2)
        for i =1:4
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
        UniquePtEdges_new = add_edge(UniqueNodes,UniquePtEdges_new);
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
        
        ind2_1 = P_P2(P_P2==P1);
        ind1_2 = P_P1
        ind1_2(ind1_2~=P2)=0;
        
        ind1_2(edge1)=0;
        ind1_2 = ind1_2(find(ind1_2))
        if size(unique(P_P1(P_P1~=0)),2)==1|P2==ind1_2
            if P2==ind1_2
                P_P2(size_p(P2,:)==max(size_p(P2,:)))=0;
                P_P2(P_P2~=P1)=0
                UniquePtEdges_new3 = cell_all{P2,1}{1,find(P_P2)};
            else
                P_P1(size_p(P1,:)==max(size_p(P1,:)))=0
                UniquePtEdges_new3 = cell_all{P1,1}{1,find(P_P1)};
                UniquePtEdges_new3 = flip(flip(UniquePtEdges_new3,1),2);
            end
            UniquePtEdges_new = [UniquePtEdges_new;UniquePtEdges_new3];
        elseif intersect(P_P1,P_P2)==0
            P1_P2 = intersect(P_P1,P_P2)
            P4 = inpoint_all(P2,find(P_P2));
            P4 = P4(1);
            UniquePtEdges_new23 = cell_all{P2,1}{1,find(P_P2)};
            P_P1_P4 = inpoint_all(P4,:);
            indP = find(P_P1_P4~=P2&P_P1_P4~=0);
            Pend = unique(P_P1_P4(indP));
            UniquePtEdges_newend = cell_all{P4,1}{1,indP};
            sizep = size_p(P1,:);
            sizep(edge1)=0;
            ind = intersect(find(P_P1==1),find(sizep==max(sizep),1));
            if size(unique(P_P1(P_P1~=0)),2)==2&P4==P1
                UniquePtEdges_new = [UniquePtEdges_new;UniquePtEdges_new23];
            else
                if ~isempty(ind)&isempty(intersect(inpoint_all(P2,:),Pend))
                    UniquePtEdges_new3 = cell_all{P1,1}{1,ind};
                    UniquePtEdges_new3 = flip(flip(UniquePtEdges_new3,1),2);
                    UniquePtEdges_new = [UniquePtEdges_new;UniquePtEdges_new23;UniquePtEdges_newend;UniquePtEdges_new3];
                else
                    UniquePtEdges_new = [UniquePtEdges_new;UniquePtEdges_new23];
                end
            end
        elseif isempty(intersect(P_P1,P_P2))
            P_P1(find(P_P1==edge1))=0;
            P_P1_P2 = [P_P1,P_P2];
            P_P1_P2(find(P_P1_P2==0))=[];
            P_P1_P2 = unique(P_P1_P2);
            P_P1_P2(P_P1_P2==P1|P_P1_P2==P2)=[];
            if size(P_P1_P2,2)==2
                P3 = P_P1_P2(1);
                P4 = P_P1_P2(2)
                if ismember(P_P1_P2(2),inpoint_all(P3,:))
                    sizep = size_p(P3,:);
                    ind32 = find(inpoint_all(P3,:)==P2)
                    if ~isempty(ind32)
                        ind23 = find(sizep==max(sizep(find(inpoint_all(P3,:)==P2))));
                        UniquePtEdges_new23 = cell_all{P3,1}{1,ind23};
                    else
                        UniquePtEdges_new23 = [];
                    end
                    
                    ind13 = find(inpoint_all(P3,:)==P1);
                    if ~isempty(ind13)
                        sizep2 = size_p(P1,:);
                        ind13 =find(sizep2==max(sizep2(ind13)));
                        UniquePtEdges_new13 = cell_all{P3,1}{1,ind13};
                    else
                        UniquePtEdges_new13=[];
                    end
                    
                    ind24 = find(inpoint_all(P4,:)==P2);
                    if ~isempty(ind24)
                        sizep2 = size_p(P4,:);
                        ind24 =find(sizep2==max(sizep2(ind24)));
                        UniquePtEdges_new24 = cell_all{P4,1}{1,ind24};
                    else
                        UniquePtEdges_new24=[];
                    end
                    
                    ind34 = find(inpoint_all(P4,:)==P3);
                    if ~isempty(ind34)
                        sizep2 = size_p(P4,:);
                        ind34 =find(sizep2==max(sizep2(ind34)));
                        UniquePtEdges_new34 = cell_all{P4,1}{1,ind34};
                    else
                        UniquePtEdges_new34=[];
                    end
                    
                    ind14 = find(inpoint_all(P4,:)==P1);
                    if ~isempty(ind14)
                        sizep2 = size_p(P4,:);
                        ind14 =find(sizep2==max(sizep2(ind14)));
                        UniquePtEdges_new14 = cell_all{P4,1}{1,ind14};
                    else
                        UniquePtEdges_new14=[];
                    end
                    
                    UniquePtEdges_new = [UniquePtEdges_new;UniquePtEdges_new13;UniquePtEdges_new23;UniquePtEdges_new24;UniquePtEdges_new34;UniquePtEdges_new14];
                end
            elseif size(P_P1_P2,2)>2
                P3 = P_P1_P2(1);
                P4 = P_P1_P2(2);
                P5 = P_P1_P2(3)
                sizep = size_p(P3,:);
                ind23 = find(sizep==max(sizep(find(inpoint_all(P3,:)==P1))));
                UniquePtEdges_new23 = cell_all{P3,1}{1,ind23};
                
                ind35 = find(inpoint_all(P5,:)==P3);
                UniquePtEdges_new34 = cell_all{P5,1}{1,ind35};
                
                ind15_egde = find(P_P2==P5);
                sizep2 = size_p(P2,:);
                ind15 =find(sizep2==max(sizep2(ind15_egde)));
                UniquePtEdges_new14 = cell_all{P2,1}{1,ind15};
                UniquePtEdges_new = [UniquePtEdges_new;UniquePtEdges_new23;UniquePtEdges_new34;UniquePtEdges_new14];
                
            end
        else
            P1_P2 = intersect(P_P1,P_P2);
            P1_P2(P1_P2==0)=[];
            P3 = P1_P2;
            P3 = P3(1);
        end
        
        if  size(unique(P_P1(P_P1~=0)),2)~=1& ~isempty(intersect(P_P1(P_P1~=0),P_P2(P_P2~=0)))
            E_P2_P3 = find(inpoint_all(P2,:)==P3);
            ind23 = find(size_p(P2,:)==max(size_p(P2,E_P2_P3)))
            UniquePtEdges_new23 = cell_all{P2,1}{1,ind23};
            E_P1_P3 = find(inpoint_all(P1,:)==P3);
            UniquePtEdges_new3 = cell_all{P1,1}{1,E_P1_P3};
            UniquePtEdges_new3 = flip(flip(UniquePtEdges_new3,1),2);
            UniquePtEdges_new = [UniquePtEdges_new;UniquePtEdges_new23;UniquePtEdges_new3];
        end
        
    end
end

return

