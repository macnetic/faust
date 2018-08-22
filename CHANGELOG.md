# FAuST Toolbox -- Flexible Approximate Multi-Layer Sparse Transform


---
**NOTE**: below this line are old changelog messages from the previous svn repository (in each start of line is displayed the svn revision and in parenthesis the corresponding git SHA1 prefix obtained after migration from svn to git).

- version 1121 (ce24328c): debut compatibility with complex dans le wrapper Matlab


- version 1116 (cfcdaa01):
	* debut compatibility with complex
	* test unitaire


- version 1115 (c6fbde76): Faust with MatGeneric
	* wrapper matlab OK
	* wrapper python Ok (but not tested on VM)


- version 1105 (ec5c6afe):
	* MATLAB wrapper (debut compatibility Faust::MatGeneric)

- version 1104 (322254ff):
	* avoid memory leak due to std::vector<Faust::MatGeneric*>, solution: create virtual destructor for abstract class Faust::MatGeneric


- version 1103 (5310d86d):
	 * test MATFILE marche


- version 1100 (ce26754e):
	* début [trunk] merge avec branch Faust_mat_generic: les faust sont une liste de matrice generique, pas forcement des matrice creuses
	* wrapper matlab HS
	* wrapper python HS

- version 1096 (90f982ed):
	* wrapper matlab OK
	* wrapper python OK
