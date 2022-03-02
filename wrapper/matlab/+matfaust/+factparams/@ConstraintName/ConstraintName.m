%% class ConstraintName
%%%

% =========================================================
%> @brief This class defines the names for the sub-types of constraints into the ConstraintGeneric hierarchy of classes.
%>
%> The table <a href="constraint.png">here</a> is a summary of the available constraints.
%>
% =========================================================
classdef ConstraintName
	properties (SetAccess = public, Constant)
		SP = 0
		SPCOL = 1
		SPLIN = 2
		NORMCOL = 3
		SPLINCOL = 4
		CONST = 5
		SP_POS = 6
		BLKDIAG = 7
		SUPP = 8
		NORMLIN = 9
		TOEPLITZ = 10
		CIRC = 11
		HANKEL = 12
		SKPERM = 13
		ID = 14
	end
	properties(SetAccess = public)
		name
	end
	methods

		function cons_name = ConstraintName(name)
			import matfaust.factparams.ConstraintName
			if(ischar(name) || iscell(name))
				name = ConstraintName.str2name_int(name);
			end
			if(name > ConstraintName.ID || name < ConstraintName.SP) %|| name == ConstraintName.BLKDIAG)
				msg = 'name must be an integer among ConstraintName.SP, ConstraintName.SPCOL, ConstraintName.NORMCOL, ConstraintName.SPLINCOL, ConstraintName.CONST, ConstraintName.SP_POS, ConstraintName.SUPP, ConstraintName.NORMLIN, ConstraintName.TOEPLITZ, ConstraintName.CIRC, ConstraintName.HANKEL, ConstraintName.SKPERM.';
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
				|| obj.name == obj.SP_POS || obj.name == obj.SPLINCOL || obj.name == obj.SKPERM;
			% BLKDIAG is a int constraint according to what faust_ConstraintGeneric.cpp indicates
		end

		function is_real = is_real_constraint(obj)
			is_real = obj.name == obj.NORMCOL || obj.name == obj.NORMLIN;
		end

		function is_mat = is_mat_constraint(obj)
			is_mat = obj.name == obj.SUPP || obj.name == obj.CONST || obj.name == obj.CIRC || obj.name == obj.TOEPLITZ || obj.name == obj.HANKEL || obj.name == obj.BLKDIAG || obj.name == obj.ID;
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
				case obj.SKPERM;
					str = 'skperm';
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
				case obj.CIRC
					str = 'circ';
				case obj.TOEPLITZ
					str = 'toeplitz';
				case obj.HANKEL
					str = 'hankel';
				case obj.BLKDIAG;
					str = 'blockdiag';
				case obj.ID
					str = 'id';
				otherwise
					error('Unknown name')
			end
		end
	end
	methods(Static)
		function id = str2name_int(str)
			import matfaust.factparams.ConstraintName
			err_msg = 'Invalid argument to designate a ConstraintName.';
			if(~ ischar(str) && ~ iscell(str))
				error(err_msg)
			end
			switch(str)
				case 'sp'
					id = ConstraintName.SP;
				case 'splin'
					id = ConstraintName.SPLIN;
				case 'spcol'
					id = ConstraintName.SPCOL;
				case 'splincol'
					id = ConstraintName.SPLINCOL;
				case 'skperm'
					id = ConstraintName.SKPERM;
				case 'sppos'
					id = ConstraintName.SP_POS;
				case 'normcol'
					id = ConstraintName.NORMCOL;
				case 'normlin'
					id = ConstraintName.NORMLIN;
				case 'supp'
					id = ConstraintName.SUPP;
				case 'const'
					id = ConstraintName.CONST;
				case 'circ'
					id = ConstraintName.CIRC;
				case 'toeplitz'
					id = ConstraintName.TOEPLITZ;
				case 'hankel'
					id = ConstraintName.HANKEL;
				case 'blockdiag'
					id = ConstraintName.BLKDIAG;
				case 'id'
					id = ConstraintName.ID;
				otherwise
					error(err_msg)
			end
		end
	end
end
