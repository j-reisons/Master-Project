function [X,numindX] = contract(X,indX,Y,indY)
% Returns contracted tensor and it's number of indices
% X is array, indX is vector containing indeces to be contracted
% Y is array, indY is vector containing indeces to be contracted
% This is the function that we stole from that other paper. It's actually
% not that hard to understand, but it's not very interesting either.

numindX = length(size(X));
numindY = length(size(Y)); % Total index counts

Xsize=size(X); % Size vectors
Ysize=size(Y);

indXr=1:numindX; 
indXr(indX)=[] ;
indYr=1:numindY ;
indYr(indY)=[]; % indXr and indYr are vectors of uncontracted indices
sizeXr=Xsize(indXr); 
sizeX=Xsize(indX); 
sizeYr=Ysize(indYr); 
sizeY=Ysize(indY); % size of contracted and remaining indices
if any(sizeX ~= sizeY) % If any contracted dimension mismatch
    error('indX and indY are not of same dimension.');
end
if isempty(indYr)
    if isempty(indXr) % If no uncontracted indeces for both
        X=permute(X,[indX]);
        X=reshape(X,[1,prod(sizeX)]);
        %Arrange them by correct contraction indeces, then reshape into
        %vector
        Y=permute(Y,[indY]);
        Y=reshape(Y,[prod(sizeY),1]);
        %Output 
        X=X*Y;
        numindX=1;
        return;
    else % If no uncontracted indeces for Y but not for X
        X=permute(X,[indXr,indX]);
        X=reshape(X,[prod(sizeXr),prod(sizeX)]);
        Y=permute(Y,[indY]);
        Y=reshape(Y,[prod(sizeY),1]);
        X=X*Y;
        Xsize=Xsize(indXr);
        X=reshape(X,[Xsize,1]);
        return
    end
end

% If both have uncontracted indeces
X=permute(X,[indXr,indX]); % Concatenated array indXl,indX. Sends summed indeces to last places
X=reshape(X,[prod(sizeXr),prod(sizeX)]); % Matrix reshape
Y=permute(Y,[indY,indYr]); % Sends summed to first places
Y=reshape(Y,[prod(sizeY),prod(sizeYr)]); % Matrix reshape
X=X*Y;
Xsize=[Xsize(indXr),Ysize(indYr)]; % Size vector of final tensor
numindX=length(Xsize); 
X=reshape(X,[Xsize,1]); % Tensors