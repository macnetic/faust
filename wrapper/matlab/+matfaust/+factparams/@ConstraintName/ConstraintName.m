%% class ConstraintName
%%%

classdef ConstraintName
	properties (SetAccess = public, Constant)
		SP = 0
		SPCOL = 1
		SPLIN = 2
		NORMCOL = 3
		SPLINCOL = 4
		CONST = 5
		SP_POS = 6
		% BLKDIAG = 7
		SUPP = 8
		NORMLIN = 9
	end
	properties(SetAccess = public)
		name
	end
	methods

		function cons_name = ConstraintName(name)
			import matfaust.factparams.ConstraintName
			if(name > ConstraintName.NORMLIN || name < ConstraintName.SP) %|| name == ConstraintName.BLKDIAG)
				msg = 'name must be an integer among ConstraintName.SP, ConstraintName.SPCOL, ConstraintName.NORMCOL, ConstraintName.SPLINCOL, ConstraintName.CONST, ConstraintName.SP_POS, ConstraintName.SUPP, ConstraintName.NORMLIN';
				error(msg)
			end
			cons_name.name = name;
		end

		function is_int = is_int_constraint(obj)
			%import matfaust.factparams.ConstraintName
			% use obj instead of ConstraintName to avoid to import the class
			% obj has access to static attributes of its class
			% (doing the same for is_real_constraint(), is_mat_constraint(), conv2str())
			is_int = obj.name == obj.SP || obj.name == obj.SPLIN || obj.name == obj.SPCOL ...
				|| obj.name == obj.SP_POS || obj.name == obj.SPLINCOL;
			% BLKDIAG is a int constraint according to what faust_ConstraintGeneric.cpp indicates
		end

		function is_real = is_real_constraint(obj)
			is_real = obj.name == obj.NORMCOL || obj.name == obj.NORMLIN;
		end

		function is_mat = is_mat_constraint(obj)
			is_mat = obj.name == obj.SUPP || obj.name == obj.CONST;
		end

		function str = conv2str (obj)
			switch obj.name
				case obj.SP
					str = 'sp';
				case obj.SPLIN;
					str = 'splin';
				case obj.SPCOL;
					str = 'spcol';
				case obj.SPLINCOL;
					str = 'splincol';
				case obj.SP_POS;
					str = 'sppos';
				case obj.NORMCOL;
					str = 'normcol';
				case obj.NORMLIN;
					str = 'normlin';
				case obj.SUPP;
					str = 'supp';
				case obj.CONST;
					str = 'const';
				%case obj.BLKDIAG;
				%	str = 'blkdiag'
				otherwise
					error('Unknown name')
			end
		end

	end
end
