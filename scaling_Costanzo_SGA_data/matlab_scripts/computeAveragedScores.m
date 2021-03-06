function averaged = computeAveragedScores(unaveraged,averagingVar)
%returns a new variable with averaged scores (over all replicates). In
%cases where only one replicate exists, a pseudo-averaging is performed.
%
%written by Sean Collins (2006) as part of the EMAP toolbox

[r c] = size(unaveraged.data);
numrep=zeros(r,1);
repID=zeros(r,10)*NaN;

% Renaming row and column labels to include the mutation type
fn=fieldnames(unaveraged);
if ismember('rowMut',fn) & ismember('colMut',fn)
    for i=1:r
        unaveraged.rowlabels(i)=cellstr([char(unaveraged.rowlabels(i)) '>>>' char(unaveraged.rowMut(i))]);
    end
    for i=1:c
        unaveraged.collabels(i)=cellstr([char(unaveraged.collabels(i)) '>>>' char(unaveraged.colMut(i))]);
    end
end

if nargin<2
    averagingVar=estimatePartnerScores(unaveraged);   %Computes the estimated average for unpaired data
    averagingVar.fitx=round(averagingVar.fitx*10);    %To make later steps work more easily
end
if isstruct(averagingVar)
    avgThresh1=min(averagingVar.fitx);
    avgThresh2=max(averagingVar.fitx);  %used to limit estimated partner scores to the range on which
                                        %they were computed for
end
% Find the number and locations of all replicate strains
orflist=union(unaveraged.rowlabels,unaveraged.collabels);
intlist=intersect(unaveraged.rowlabels,unaveraged.collabels);
numTot=length(orflist);
natPresent=zeros(size(unaveraged.rowlabels(:)));
for i=1:r
    if countnan(unaveraged.data(i,:))<c
        natPresent(i)=1;
    end
end
kanPresent=zeros(size(unaveraged.collabels(:)));
for i=1:c
    if countnan(unaveraged.data(:,i))<r
        kanPresent(i)=1;
    end
end
for i=1:numTot
    ind1=find(strcmp(orflist(i),unaveraged.rowlabels(:))& natPresent>0);
    ind2=find(strcmp(orflist(i),unaveraged.collabels(:))& kanPresent>0);
    numrep(i)=length(ind1)+length(ind2);
    repID(i,1:numrep(i))=[ind1' -1*ind2'];
end

averaged.data=zeros(numTot,numTot)*NaN;
%Now computing averages for the strains we have at least one replicate of
for i=1:numTot
    for j=1:(i-1)
        if (numrep(i)>0)&&(numrep(j)>0)
            list = getDataPoints(unaveraged,numrep,repID,i,j);
            if length(list)==1 && isstruct(averagingVar)
                score1=min(max(list(1)*10,avgThresh1),avgThresh2);
                score1=(round(score1)); 
                estPartner=averagingVar.fity(averagingVar.fitx==score1);
                averaged.data(i,j)=(estPartner + list(1))/2;   %psuedo-averaging
            elseif ~isempty(list)
                averaged.data(i,j)=mean(list);
            end
            counters(i,j)=length(list);
            averaged.data(j,i)=averaged.data(i,j);
        end
    end
end
averaged.collabels=orflist(:);
averaged.rowlabels=orflist(:);
if ismember('geneToOrf',fn)
    averaged.geneToOrf=unaveraged.geneToOrf;
end

% Compact the matrix by removing rows (and the corresponding columns) that
% have been collapsed into other averages
averaged=pruneScoreMatrix(averaged);

if length(intlist)==0   %the set of original rows and columns were disjoint
    temp=averaged;
    [averaged.rowlabels,ind1a,ind1b]=intersect(averaged.rowlabels,unaveraged.rowlabels);
    [averaged.collabels,ind2a,ind2b]=intersect(averaged.collabels,unaveraged.collabels);
    averaged.data=temp.data(ind1a,ind2a);
end

if ismember('rowMut',fn) & ismember('colMut',fn)
    for i=1:length(averaged.rowlabels)
        temp=parseString(averaged.rowlabels(i),'>>>');
        averaged.rowlabels(i)=temp(1);
        averaged.rowMut(i,1)=temp(2);
    end
    for i=1:length(averaged.collabels)
        temp=parseString(averaged.collabels(i),'>>>');
        averaged.collabels(i)=temp(1);
        averaged.colMut(i,1)=temp(2);
    end
end

end


%*************** EXTRA FUNCTIONS *************************

function list=getDataPoints(scoremat,numrep,repID,i,j)
    list=[];
    for i1=1:numrep(i)
        for j1=1:numrep(j)
            if (repID(i,i1)>0)&(repID(j,j1)<0)
                list=[list scoremat.data(repID(i,i1),-1*repID(j,j1))];
            end
            if (repID(i,i1)<0)&(repID(j,j1)>0)
                list=[list scoremat.data(repID(j,j1),-1*repID(i,i1))];
            end
        end
    end
    list=list(~isnan(list));
end
