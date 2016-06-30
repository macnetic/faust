function [Ddvp, cum_dvp] = cum_dvp(D)
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
cum_dvp=cell(1,Nlayers);

if Nlayers==0
    Ddvp = 1;
else    
    cum_dvp{1} =D{1};
    for i = 2:Nlayers
        if size(D{i},1) == size(cum_dvp{i-1},2)
            cum_dvp{i} = cum_dvp{i-1} * D{i};
        end
    end

    Ddvp=cum_dvp{1,Nlayers};
end

end

