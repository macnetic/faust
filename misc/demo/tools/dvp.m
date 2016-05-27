function Ddvp = dvp(D)
%dvp. Development of the input factorized matrix.
%  Ddvp = dvp(D) develops the cell-array of matrices D into the matrix Ddvp
%  which is the product of the matrices contained in D:
%  Ddvp = D{1}*D{2}*...*D{n}.
%
%  References:
%  [1] Le Magoarou L. and Gribonval R., "Learning computationally efficient
%  dictionaries and their implementation as fast transforms", submitted to
%  NIPS 2014

%  Luc Le Magoarou
%  PANAMA team
%  Inria, Rennes - Bretagne Atlantique
%  luc.le-magoarou@inria.fr
%
%  June 2014

Nlayers = size(D,2);
if Nlayers==0
    Ddvp = 1;
else    
    Ddvp =D{1};
    for i = 2:Nlayers
        if size(D{i},1) == size(Ddvp,2)
            Ddvp = Ddvp * D{i};
        end
    end
end
end