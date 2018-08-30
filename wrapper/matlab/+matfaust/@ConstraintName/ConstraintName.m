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
			import matfaust.ConstraintName
			if(name > ConstraintName.NORMLIN || name < ConstraintName.SP) %|| name == ConstraintName.BLKDIAG)
				msg = 'name must be an integer among ConstraintName.SP, ConstraintName.SPCOL, ConstraintName.NORMCOL, ConstraintName.SPLINCOL, ConstraintName.CONST, ConstraintName.SP_POS, ConstraintName.SUPP, ConstraintName.NORMLIN';
				error(msg)
			end
			cons_name.name = name;
		end

		function is_int = is_int_constraint(obj)
			import matfaust.ConstraintName
			is_int = obj.name == ConstraintName.SP || obj.name == ConstraintName.SPLIN || obj.name == ConstraintName.SPCOL ...
				|| obj.name == ConstraintName.SP_POS || obj.name == ConstraintName.SPLINCOL;
			% BLKDIAG is a int constraint according to what faust_ConstraintGeneric.cpp indicates
		end

		function is_real = is_real_constraint(obj)
			import matfaust.ConstraintName
			is_real = obj.name == ConstraintName.NORMCOL || obj.name == ConstraintName.NORMLIN;
		end

		function is_mat = is_mat_constraint(obj)
			import matfaust.ConstraintName
			is_mat = obj.name == ConstraintName.SUPP || obj.name == ConstraintName.CONST;
		end

		function str = conv2str (obj)
			import matfaust.ConstraintName
			switch obj.name
				case ConstraintName.SP
					str = 'sp';
				case ConstraintName.SPLIN;
					str = 'splin';
				case ConstraintName.SPCOL;
					str = 'spcol';
				case ConstraintName.SPLINCOL;
					str = 'splincol';
				case ConstraintName.SP_POS;
					str = 'sppos';
				case ConstraintName.NORMCOL;
					str = 'normcol';
				case ConstraintName.NORMLIN;
					str = 'normlin';
				case ConstraintName.SUPP;
					str = 'supp';
				case ConstraintName.CONST;
					str = 'const';
				%case ConstraintName.BLKDIAG;
				%	str = 'blkdiag'
				otherwise
					error('Unknown name')
			end
		end

	end
end
