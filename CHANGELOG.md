# FAuST Toolbox -- Flexible Approximate Multi-Layer Sparse Transform

(for further details take a look to full git logs)

- f0858790 Update gen_matlab_inline_doc_from_doxy_blocks.sh to make it work on macOS.
- 60cc5a78 **Update API doc (py/mat).**
- 22b4f51a Change display function for a Faust.
- c1d5ac23 **Enhance random Faust generation.**
- c0ae3a3f Disable gen_matlab_inline_doc_from_doxy_blocks.sh on macOS because currently the script works only on Linux.
- 182aa101 **Add expire dates for packages sent to gitlab and keep a local copy of revision pkg for macOS.**
- 10092f2e **Add equivalent functions matfaust.Faust.numel() and pyfaust.Faust.size().**
- 565c57bc Fix in cmake script a bug introduced by e92d4bc6.
- 9dbfffea Fix missing namespace matfaust in matfaust.Faust.mtimes_trans() (next to 3f4d86cc).
- ef7efd74 Update post-install script for py. wrapper (in next to edfcad65).
- 06962bd0 Check that bounds are respected in matfaust.Faust.subsref().
- 3f4d86cc **Create a namespace for matlab wrapper: matfaust.**
- e92d4bc6 Add a script for editing Faust.m to generate a Matlab inline doc from doxygen code doc.
- 14aca00e **(tag: 2.2rc20) Bind the C++ faust mul. (cb58f3a9) into matlab wrapper Faust.mtimes().**
- cb58f3a9 **Implement the mul. of two Fausts in C++ and bind py. wrapper to this impl.**
- 1c661768 **(tag: 2.2rc19) Bind matlab Faust.susbref() to the new impl. of slicing/subscript (c4c2c431).**
- 96ca57d9 Fix Faust::MatSparse::get_rows/cols(): nnz wasn't updated in submatrix.
- c4c2c431 **(tag: 2.2rc18) Reimplement completely pyfaust.Faust.__getitem__(): now it returns a Faust object instead of a ndarray.**
- 918b4a41 (tag: 2.2rc17) Update API doc.
- edfcad65 **Rename module FaustPy to pyfaust.**
- 1817ea89 Replace FaustPy.Faust.isReal() by a dtype property.
- 2fe19fce **Add property decorated functions FaustPy.Faust.H() and FaustPy.Faust.T().**
- 8b9969b2 Add postinst. script for macOS to open the readme automatically after installation.
- d9eafcb0 Set Faust.rand() more flexible in both wrappers by changing arguments (polymorphism).
- 821b2cd8 Temporarly disable Linux nightly tests.
- cf974102 Rename matlab wrapper Faust.get_fact() to Faust.get_factor().
- 7bb8e123 Exclude Matlab wrapper functions mtimes_trans and subsasgn from doxygen doc.
- 1d684e27 Rename Faust.todense() to Faust.toarray() and add a new Faust.todense() which returns a numpy.matrix.
- c0749139 Update API doc.
- b10d1da1 (tag: 2.2rc16) Rename FaustPy.Faust.get_nb_factors() to get_num_factors() and likewise for matlab wrapper get_nb_factor().
- c4467e97 Modify the way identity matrix is initialized in cpp unit tests for mul.
- 23385153 Fix wrong display of Faust transpose and update API doc.
- 7b68675f Override default \_\_repr\_\_() in Faust class and use it to refactor display().
- 8e449471 Add a new function to_string() in Faust::MatGeneric/MatSparse/MatDense/Transform to refactor their Display() functions and be able to get information string from a Faust without having to display it.
- 13e61646 Delete FaustPy.Faust.get_nb_rows()/get_nb_cols().
- 1baa4826 Rename FaustPy.Faust.size() to shape and use a property (decorator) instead of function.
- 6f1e7d5a Rename nnz() to nnz_sum() (in matlab and py wrappers).
- c8c6b182 Minor fix in test script.
- 476878a0 (tag: 2.2rc15) Rename RCG() to rcg() (matlab and py wrappers).
- fa0b6233 **(tag: 2.2rc14) Update guide doc.**
- aa0d2d05 (tag: 2.2rc13) Add a new cmake option (and impl.) to generate a README with install. instructions for all systems.
- 86a85697 Install macOS release version into prefix /opt/local/faust like for linux (removing tag from folder name).
- ca70e3f8 **(tag: 2.2rc12) Update readmes and their generation cmake script.**
- 374abc73 (tag: 2.2rc11) Limit the gitlab-pages directory to doc/html (not the whole doc anymore).
- 5930d9af Update doxygen mainpage.md
- e8a9ac21 Update macOS installer readme.
- 3abc3e07 **(tag: 2.2rc10) Send packages linux and mac to gitlab.**
- 821dec12 Update API doc.
- 575abeb7 Add an option to pick operator norm or Frobenius to compute relative error. Calculate the norm directly with singular values instead of computing the approximate of M and then calculate the norm with it.
- d7271db4 Update API/code doc.
- e315a460 Plot RC/density instead of RCG.
- 23cf8c94 (tag: 2.2rc9) Fix bug in Faust::Transform<FPP,Cpu>::updateNonZeros().
- c251430c Update API/code doc.
- e7ab567a Fix python unit tests for RCG and density next to the fix c586c85b.
- 98d37886 Fix constructor ignoring of multiplicative scalar when Faust obj. is constructed from a file (matlab wrapper).
- e1b20d5b Add a faust test script showing how the .mat file size is related to RCG.
- c586c85b Fix Faust.density() and RCG() for Py2.
- 46d60f4b **Implement multiplicative constant to apply to a Faust at initialization time.**
- 9fc63bdb **Handle also sparse matrices (csr_matrix) in list_factors argument of constructor Faust().**
- f65784dd Fix wrong index bug in one constructor Faust::MatSparse().
- e05964fd Add a script for showing truncated SVD accuracy versus the RCG.
- 79463a12 **Add gitlab pages with doxygen doc for each tag.**
- 28e87505 Add pkg postinst. script to symlink matio if needed on debian/ubuntu sys.
- 86abf82d Shorten the commit sha to its unique eight first digits in package names and versions.
- 4c5fe261 Replace in doxygen doc the static date by the date of generation and display the cpack package version.
- f46b6792 Enhance deb/rpm post-installation script for matlab wrapper (more info and verifs).
- c01b9af9 Enhance deb/rpm post-installation script for python wrapper (more info and verifs).
- c2005f87 Generalize matio static libs option to matlab in addition to python (143be5e4).
- b338c711 **Add a new package generation for linux using static libs for python wrapper.**
- e30641ba Refactor raw coded data about matio from setup.py.in to CMakeLists.txt.
- 7b2fb435 Compile and test for python2 even if python3 is found and likewise search the two versions of python and cython.
- e8181d22 Add missing description for RPM package (based on debian pkg description).
- 4727b49d Replace URL http://faust.gforge.inria.fr by http://faust.inria.fr in a lot of files.
- de00d5b9 Avoid to exit in matlab postinstall script.
- 4e120fbe Add a new postinstall script specific to macOS.
- c601e7b6 Remove /opt from auto filelist for rpm package because it's a conflict with filesystem package.
- 5f3358c3 Specify dependencies for rpm package and filter matlab's.
- aeada65a Update py wrapper doc.
- cd18a6a2 Filter more properly the doxygen doc when excluding C++ doc/static lib to keep only Python and/or Matlab wrapper doc.
- 5c260be6 Add package postinstall script to automatize py wrapper setup.
- f98e2eea Fix memory leak in Faust::MatSparse::randMat().
- d328d728 **Add Faust.rand() to the matlab wrapper.**
- fe52e009 Update py. wrapper code doc.
- 2463997b Fix two randFaust() bugs on MacOSX/clang/py2.7.
- 62b51ecd Fix bug in py wrapper display() occurring for complex Fausts.
- c99d599a **Add a random faust generation method in cpp lib and py. wrapper.**
- fc80786c Enable cmake code for C++ version detection and corresponding compiler flag setup.
- f7cedb93 **Implement Frobenius norm in cpp core lib. and bind py and matlab wrappers to it. Add related unit tests.**
- 1d2e57ee **Implement sparse matio/matlab file output for Faust factors with support of real and complex scalars, Faust's conjugate and transpose.**
- f920ad4b Fix bug. Prevent save of empty Faust.
- a759047c **Add L1-norm to core cpp lib, py. and matlab wrappers. Also add unit tests for this norm for py. and matlab.**
- c2a191e6 Fix matlab Faust destructor bug occurring when factors passed are invalid because of not matching dimensions and hence the Faust object is not really instantiated.
- c748449f **Implement Faust over complex field into Python wrapper and add related unit tests.**
- 0e8c73f4 **Add Matlab wrapper implementations of ctranspose() and conj() with associated unit tests. For Python wrapper add conj() and getH() (for conjugate transpose)** but the unit tests are not relevant since we haven't yet the complex matrix support.
- 1d42602c Add unit tests for Faust deletion.
- c775e359 **Move to a weak reference design pattern for Transform objects used in wrappers. Fix memory leaks for Transform.**
- 2c75ad90 Add unit tests for Matlab wrapper Faust class.
- a68a4c6a Slightly modify unit tests for norm(), \__getitem_\_() and test params (num. of factors and dimension sizes).
- 569be809 **Handle save of complex or complex-transpose faust in Faust::MatDense::toMatIOVar().**
- 807d7845 **Handle save of a transpose Faust in Faust::Transform::save_mat_file()** and Faust::MatDense::toMatIOVar()/Faust::MatSparse::toMatIOVar(). Allow py and matlab wrappers to access the transpose case by passing a boolean to core lib.
- 3ab67897 **Add unit tests for not covered yet FaustPy.Faust functions.**
- d463e3b3 Remove flag CMAKE_BUILD_TYPE=Release because it fails package generation for Linux.
- 9ab90811 Add matlab doxydoc if matlab wrapper is built.
- e0e055fb **Add post-installation scripts into packages for auto-setting matlab wrapper path into a matlab startup script.**
- 49418cf5 Remove wrapper language specific code to save a Faust (Py and matlab).
- bd2edf42 **Add doxygen code doc generation for the Matlab wrapper.**
- b3f59bbc **Add a Matlab unit tests class for Faust class defined in Faust.m.**
- a5cc8a69 Fix and factorize Faust.delete().
- e0825fe2 Include licenses and README in MacOSX packages and doc directory at installation stage.
- cf95ce6f Set the linux package version accordingly to the git commit/tag.
- 151a1430 **Add a unit tests script for Python wrapper main class FaustPy.Faust.**
- 15b2d5d2 Fix bug into Faust::Transform::save_mat_file(), the number of faust factors saved was always two.
- 57d19063 Document the code by making profit of doxypypy and do other minor changes.
- 59db7796 Temporarly disable Windows nightly tests.
- da9257c9 Check that FaustCorePy.get_fact() index argument is in correct range.
- 4a827930 **Add API generated documentation files to packages.**
- 7668557b **Bind matlab wrapper to Faust::Transform::save_mat_file().**
- 4e9559dc Add libmatio dyn. lib. to MacOSX in ci package generation.
- 2aac8ad7 **Add in core cpp lib the capabililty to write a FAuST into a Matlab file and link the Py wrapper to it.**
- fc80aa72 Remove empty definitions of multiply() in MatGeneric and replace by virtual pure function members.
- 955eb586 Make also a package with pkgbuild for MacOSX.
- e78e0e16 **Add FaustPy doxygen doc. generation.**
- 3d4027ed Fix Faust.get_factor() for the transpose case and add related tests into quickstart.py.
- c56e13a3 **Add functions in the Py wrapper: size(), get_factor(), save() and adapt the constructor to read a Faust from file.**
- b44bb4c8 Update README.
- 8930b6d5 Add ftp command to upload package for linux.
- 387182db Convert README from txt to Markdown.
- d62e1945 **Add methods to Py wrapper: nnz(), density(), RCG(), norm() and get_nb_factors().** These methods are already in Matlab wrapper.
- e3628a52 Install doxygen doc for core c++ API if core library's installed.
- 94f5b2aa Strip binaries in Release build mode and optionally exclude core faust lib.
- 1ee8d7f3 Exclude git tags from non-related ci jobs.
- f301a09b Reindent the code.
- 3210b598 Correct error in R2016a/gcc versions requirements.
- a517d067 **Add ci jobs for packaging when git tagging.**
- 3d64cf65 Exclude clang/llvm from gcc/matlab versions verification.
- ae3c6ff1 Add tags and except sections to ci job for linux packaging.
- 2178f84e Add a ci job for packaging faust on Linux (only rpm and deb).
- 79e67df2 Check gcc version matches with Matlab's, block compilation if not.
- fa6c4538 Update MacOS packaging job.
- 18aee7aa Exclude scheduled ci pipelines from revision packaging job and add the commit sha as a suffix to package name.
- 2bd60833 Add switch to sudo in order to get pass from stdin.
- 61fbf78b Add a new ci job for packaging faust on MacOS X.
- e7834f85 Avoid plotting faust time for py wrapper test when we are in test env. (with no display).
- 489f5011 Fix bug with a matrix-vector product into py wrapper test.
- c2fc2528 Fix 2 minor warnings in cmake scripts.
- 2e42b46a Fix numpy array access bug, issue #9.
- f4ddfee1 Add gitlab ci job for Win7.
- 3c227222 Update CDashConfScript.cmake for Windows.
- 0a0f191e Decrease ci test time and add specific jobs for macos and linux nightly tests.
- 6e8aa606 Change installation prefix directory for testing.
- 1ba46955 Update CDashConfScript.
- 10864594 Remove useless command ctest_empty_binary_directory().
- 9d994256 Remove useless ci job tags.
- 25d47f9a Rename CDashConfScriptLinux.cmake to CDashConfScript.cmake.
- bfec94a6 Separate gitlab ci jobs into quick and slow jobs and add job tags for win and macos.
- e64d0c66 Fix clang compilation error in matlab wrapper.
- 7b749138 Fix test related cmake error on Windows.
- f64202aa Fix matlab mex compilation error on Windows.
- c87d75ff Fix bug with matlab environment probing on Windows.
- 5f4092c7 Add python and matlab wrappers to linux tests.
- 5134963c Update in order to make python3 available.
- 41312154 Fix gcc/linux compilation errors for types uint32/64_t declared twice.
- 46586264 Add CDash/CTest site name.
- de5fffe4 Add ci yml script and ctest/cdash script called by it.

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
