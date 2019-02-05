function scoremat = pruneScoreMatrix(scoremat)
%removes rows and columns which contained only NaN values in the data
%matrix
%
%written by Sean Collins (2006) as part of the EMAP toolbox


[r c]=size(scoremat.data);
fn=fieldnames(scoremat);

ind1=[];
ind2=[];
for i=1:r
    if countnan(scoremat.data(i,:))==c
        ind1=[ind1 i];
    end
end
for i=1:c
    if countnan(scoremat.data(:,i))==r
        ind2=[ind2 i];
    end
end

scoremat.rowlabels(ind1)=[];
scoremat.collabels(ind2)=[];
scoremat.data(ind1,:)=[];
scoremat.data(:,ind2)=[];
if ismember('rowCoord',fn)
    scoremat.rowCoord.chrom(ind1)=[];
    scoremat.rowCoord.coord(ind1,:)=[];
    scoremat.colCoord.chrom(ind2)=[];
    scoremat.colCoord.coord(ind2)=[];
end
if ismember('rowMut',fn)
    scoremat.rowMut(ind1)=[];
    scoremat.colMut(ind2)=[];
end
