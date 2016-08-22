function A_normalize = normalizeCol(A,varargin)

nb_input=length(varargin);

if (nb_input == 1)
    lambda=varargin{1};
else
    lambda=1;
end

if(iscell(A))
    X =lambda*dvp(A);
else
    X = A;
end

X(:,sum(X.^2,1)==0)=1;% null column are fixed to 1

if (iscell(A))
    A_normalize=A;
    A_normalize{end}=A_normalize{end}./repmat(sqrt(sum(X.^2,1)),size(X,1),1);
else
    A_normalize = A./repmat(sqrt(sum(X.^2,1)),size(X,1),1);
end

end

