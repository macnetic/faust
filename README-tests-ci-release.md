# Automatic testing and continuous integration of FAµST with Gitlab

This document briefly describes the tests, continuous integration
and maintenance related elements for the FAµST project.
There are specific tests for each package/language used in the project. You'll
find a section for each one of them in the next of this document.

### Summary
1. [Test policy on Gitlab](#gitlab_pol_test)  
2. [Testing pyfaust/matfaust](#py_matfaust_tests)  
2.1 [Test scripts for pyfaust/matfaust](#py_matfaust_tests)  
2.2 [Gitlab continuous integration (ci) job](#py_mat_test_ci_jobs)  
2.3 [Test report on Gitlab pages](#py_mat_test_reports)  
2.4 [Doctest](#doctest)  
2.5 [Test coverage](#test_cover)
3. [FAµST C++ core tests](#cpp_tests)  
3.1 [Tests](#cpp_tests)  
3.2 [Gitlab continuous integration (ci) job](#py_mat_test_ci_jobs)  
3.3 [Test report on Gitlab pages](#cpp_test_report)
4. [Package tests and automatic release](#test_release_packages)  

<a name="gitlab_pol_test"/>

**1. Test policy on Gitlab**

A strongly strict policy is applied on Gitlab test jobs.
The goal is to ensure the non-regression of the project (which has a large base of
code) and limit avoidable bugged release.
In practice, if any test, unit test or doctest fails during a test job then it won't pass
on Gitlab, the pipeline will be marked as a failed pipeline (so it is reliable
and instantaneous to verify all's ok).
The failed tests are clearly indicated in the job output or into the test report
(see below for how to find the test reports on gitlab).
A Gitlab pipeline failure means no package will be generated or at least that
they won't/shouldn't reach the release stage.

**Note**: it is possible for a maintainer of the project to ignore some tests
results and push a git-tag to release new packages but in that case at least
the packages tests must pass or it will fail).

**2. Testing pyfaust/matfaust**  


**<a name="py_matfaust_tests">
2.1 Test scripts for pyfaust/matfaust
</a>**

- Directory location: ``misc/test/src/Python`` for pyfaust and
  ``misc/test/src/Matlab``
  for matfaust.

**Note**: at the time of writing, only 5/18 of the pyfaust's test scripts and 8/13
of matfaust's are integrated in the automatic tests. The other scripts cover other
needs (e.g. benchmark, debug, etc.).
Because these scripts are old there is no guaranty they still work (maybe an update to
recent FAµST API is needed).

- Unit tests:
    * ``test_FaustPy.py`` (124 unit tests at the time of writing for
      ``pyfaust.Faust`` class and all algorithms/functions),
    * ``FaustTest.m`` (31 unit tests for the ``matfaust.Faust`` class), ``FaustFactoryTest.m`` (39 unit tests for factorizations and other algorithms).

**<a name="py_mat_test_ci_jobs">
2.2 Gitlab continuous integration (ci) job
</a>**

- All ci jobs are defined in the ``.gitlab-ci.yml`` file.
- Gitlab ci jobs: ``ctest_python`` for pyfaust, ``ctest_matlab`` for matfaust.
- how it works:
  [ctest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html)
  is triggered by the ci job and after a building step it
  executes (on a Linux/MacOS VM) the test and unit tests above and other test scripts.
  The configuration of ctest for running the tests is located in ``misc/test/CMakeLists.txt``

**<a name="py_mat_test_reports">
2.3 Test reports on Gitlab pages
</a>**

- On the end of ``ctest_python``/``ctest_matlab`` execution a test report in HTML
  is automatically uploaded to the project Gitlab pages (as a hidden page).
  The link is displayed in the ci job output (the finished jobs are listed
  [here](https://gitlab.inria.fr/faustgrp/faust/-/jobs)).
  Here is [an example of report for pyfaust](https://faustgrp.gitlabpages.inria.fr/-/faust/-/jobs/3163959/artifacts/build_FaustLinuxPython/python_pyfaust_test_output.html)
  (if this link doesn't work, the associated pipeline has been deleted on Gitlab,
  please look [here](https://gitlab.inria.fr/faustgrp/faust/-/jobs) for a recent ``ctest_python``/``ctest_matlab`` output).

**<a name="doctest">
2.4 Doctest
</a>**

The pyfaust/matfaust API comes with many examples integrated on the [documentation](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/index.html).
They are all fully tested in the
``doctest_nightly_pyfaust``/``doctest_nightly_matfaust`` ci job (cf. ``gitlab-ci.yml``)
and comply with the job passing policy exposed in [1.](#gitlab_pol_test) (any doctest failure means the ci pipeline
fails). These ci jobs are [scheduled](https://gitlab.inria.fr/faustgrp/faust/-/pipeline_schedules)
to run every night.
If any of the tested pyfaust/matfaust submodules or functions fails (and they are all tested at the time of
writing), the erroneous tests will be displayed in the ci job output.

**<a name="test_cover">
2.5 Test Coverage
</a>**

It really matters to keep a metric of how much code is covered by the tests in the project.

- The ci job ``pyfaust_test_code_coverage`` measures the test coverage for
pyfaust (using the same tests as ``ctest_python`` and
``doctest_nightly_pyfaust``). At the moment I'm writing these lines, the cover
is about 72%. The same ci job produces a report available on gitlab-pages
too as job's artifact (e.g.: [report](https://faustgrp.gitlabpages.inria.fr/-/faust/-/jobs/3169837/artifacts/htmlcov/index.html) -- note: this link won't last forever). The report is also available directly as text in the ci job output.

- For matfaust the ci job is ``matfaust_test_code_coverage`` it uses MOCov to
produce a coverage report as this [one](https://faustgrp.gitlabpages.inria.fr/-/faust/-/jobs/3186956/artifacts/coverage_html/index.html)
(this link might be deleted in the future, in which case you might look at artifacts of recent ci
[jobs](https://gitlab.inria.fr/faustgrp/faust/-/jobs)). The tests used for coverage calculation
are pretty much the same as the ones in ``ctest_matlab`` and ``doctest_nightly_matfaust``.
When writing this doc, the coverage was about 48.6%, it is partially due to the fact that the
code for mex float (single in matlab) is not auto-tested but there is no reason to think it
would fail more often than the code for double which is tested. However it would be advisable
to pursue for a larger coverage.

- About GPU code: unfortunately no automatic testing is made for GPU code
  because no NVIDIA GPU is available on VMs/Docker daemons that run the project gitlab-runners.

**<a name="cpp_tests">
3. FAµST C++ core tests
</a>**

### 3.1 Tests

- Directory location: ``misc/test/src/C++`` and ``misc/test/src/C++/unit``

### 3.2 Gitlab continuous integration (ci) job
<a name="cpp_test_ci_jobs"/>

- Gitlab ci job: ``ctest`` performs C++ tests.
- For more details about ctest, take a look at [2.2](#py_mat_test_ci_jobs).

### 3.3 Test report on Gitlab pages
<a name="cpp_test_report"/>

- On the end of ``ctest`` ci job execution a test report in HTML
  is automatically uploaded to the project Gitlab pages (as a hidden page).
  The link is displayed in the ci job output (the finished jobs are listed
  [here](https://gitlab.inria.fr/faustgrp/faust/-/jobs)).
  Here is [an example of report for C++ tests](https://faustgrp.gitlabpages.inria.fr/-/faust/-/jobs/3193880/artifacts/build_FaustLinux/cpp_test_report.html)
  (if this link doesn't work, the associated pipeline has been deleted on Gitlab,
  please look [here](https://gitlab.inria.fr/faustgrp/faust/-/jobs) for a recent ``ctest`` output).

**<a name="test_release_packages">
4. Package tests and automatic release
</a>**

TODO
