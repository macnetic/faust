# Automatic testing and continuous integration of FAµST with Gitlab

This document briefly describes the tests, continuous integration
and maintenance related elements for the FAµST project.
There are specific tests for each package/language used in the project. You'll
find a section for each one of them in the next of this document.

[[_TOC_]]


### Test policy on Gitlab

A strongly strict policy is applied on Gitlab test jobs.  The goal is to ensure the
non-regression of the project (which has a large code base) and limit avoidable
bugged releases.  In practice, if any test, unit test or doctest fails during a test
job then it won't pass on Gitlab, the pipeline will be marked as a failed pipeline
(so it is reliable and instantaneous to verify all's ok).  The failed tests are clearly
indicated in the job output or into the test report (see below for how to find the
test reports on Gitlab).  A Gitlab pipeline failure means no package will be generated
or at least that they won't/shouldn't reach the release stage. Below is the
status (the badge) for the last ran pipeline.

[![pipeline status](https://gitlab.inria.fr/faustgrp/faust/badges/hakim_branch/pipeline.svg?ignore_skipped=true)](https://gitlab.inria.fr/faustgrp/faust/commits/hakim_branch)

**Note**: it is possible for a maintainer of the project to ignore some tests
results and push a git-tag to release new packages but in that case at least
the packages tests must pass or it will fail.

### Testing pyfaust/matfaust/C++ core


#### Test scripts/programs

- Directory location: ``misc/test/src/Python`` for pyfaust and
  ``misc/test/src/Matlab`` for matfaust. For C++, look at ``misc/test/src/C++``.

**Note**: at the time of writing, only 5/18 of the pyfaust's test scripts and 8/13
of matfaust's are integrated in the automatic tests. The other scripts cover other
needs (e.g. benchmark, debug, etc.).
Because these scripts are old there is no guaranty they still work (maybe an update to
recent FAµST API is needed).

- Unit tests:
    * ``test_FaustPy.py`` (128 unit tests at the time of writing for
      ``pyfaust.Faust`` class and all algorithms/functions),
    * ``FaustTest.m`` (34 unit tests for the ``matfaust.Faust`` class), ``FaustFactoryTest.m``
      (39 unit tests for factorizations and other algorithms).
    *  For C++ code ``misc/test/src/C++/unit`` (but note that it's a bit misnamed
       because unit tests are also available in the parent directory). Up to 167 tests
       are available (accounting for demultiplied tests based on different scalar
       types available -- float, double, complex float, complex double).

#### Gitlab continuous integration (CI) jobs

- All CI jobs are defined in the ``.gitlab-ci.yml`` file.
- Gitlab CI jobs: ``ctest_python`` for pyfaust, ``ctest_matlab`` for matfaust. ``ctest``
  CI job performs C++ tests.
- How it works:
  [ctest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html)
  is triggered by the CI job and after a building step it
  executes (on a Linux/MacOS VM) all the unit tests above and other test scripts.
  The configuration of ctest for running the tests is located in ``misc/test/CMakeLists.txt``

-  Note that for a matter of pipeline speed the CI job named ``ctest_nightly`` is
  made to run the tests that are very slow/long.  They are enabled through the cmake
  option ``SLOW_TESTS`` (which is defaultly set to ``OFF``). There is one job version
  for each supported OS (``ctest_nightly_linux``, ``ctest_nightly_macos``, ``ctest_nightly_win10``).
  On Linux and Mac OS all the tests are performed (pyfaust, matfaust, C++ core) while
  on Windows on C++ tests are launched.  These tests are scheduled on Gitlab to run
  periodically by night. Note that also test coverage CI jobs run at night, it allows
  to ensure that we keep track of the coverage rates for pyfaust, matfaust and C++
  core (see [Test Coverage](#test-coverage).

#### Test reports on Gitlab pages

- On the end of ``ctest_python``/``ctest_matlab`` execution a test report in HTML
  is automatically uploaded to the project Gitlab pages.
  The link is displayed in the CI job output (the finished jobs are listed
  [here](https://gitlab.inria.fr/faustgrp/faust/-/jobs)).
  Here is [an example of report for pyfaust](https://faustgrp.gitlabpages.inria.fr/-/faust/-/jobs/3163959/artifacts/build_FaustLinuxPython/python_pyfaust_test_output.html)
  (if this link doesn't work, the associated pipeline has been deleted on Gitlab,
  please look [here](https://gitlab.inria.fr/faustgrp/faust/-/jobs) for a recent ``ctest_python``/``ctest_matlab``
  output). You can also simply browse the job directory of artifacts to find the report.

- On the end of ``ctest`` CI job execution, a test report in HTML is produced for
  the C++ tests too.  It is automatically uploaded to the project Gitlab pages.
  The link is displayed in the CI job output (the finished jobs are listed
  [here](https://gitlab.inria.fr/faustgrp/faust/-/jobs)).  Here is
  [an example of report for C++ tests](https://faustgrp.gitlabpages.inria.fr/-/faust/-/jobs/3193880/artifacts/build_FaustLinux/cpp_test_report.html)
  (if this link doesn't work, the associated pipeline has been deleted on Gitlab, please
  look [here](https://gitlab.inria.fr/faustgrp/faust/-/jobs) for a recent ``ctest``
  output).  At the time of writing a group of 167 C++ tests are ran by the ctest CI
  job. It includes many unit tests and all these tests passed the latest pipeline.
  You can also simply browse the job directory of artifacts to find the report.

- Alternative test reports are available here: [FAµST CDash](https://cdash-ci.inria.fr/index.php?project=faust),
 they are uploaded by ctest/cdash (CDash configuration is available in ``CDashConfScript.cmake``
 at the root of the project). It gives additional information about the
 building part of tests (which is important in case of errors, that's where you might look up
 to understand why a C++ test or a wrapper failed to build).

#### Doctest

The pyfaust/matfaust API comes with many examples integrated in the
[documentation](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/index.html).
They are all fully tested in the ``doctest_pyfaust``/``doctest_matfaust``
ci jobs (cf. ``gitlab-ci.yml``) and comply with the job passing policy exposed in [1.](#gitlab_pol_test)
(any doctest failure means normally that the CI pipeline fails).
If any of the tested pyfaust/matfaust submodules or functions fails (and they are
all tested at the time of writing), the erroneous tests will be displayed in the CI
job output.

#### Test Coverage

It really matters to keep a metric of how much code is covered by the tests in the project.

- The CI job ``pyfaust_test_code_coverage`` measures the test coverage for
  pyfaust (using the same tests as ``ctest_python`` and ``doctest_pyfaust``).
  Look at the badge below to know what is the current coverage of pyfaust tests. The
  same CI job produces a report available on gitlab-pages too as job's artifact
  (e.g.: [report](https://faustgrp.gitlabpages.inria.fr/-/faust/-/jobs/3169837/artifacts/htmlcov/index.html)
  -- note: this link won't last forever). The report is also available directly as
  text in the CI job output.

![pyfaust test coverage](https://gitlab.inria.fr/faustgrp/faust/badges/hakim_branch/coverage.svg?job=pyfaust_test_code_coverage&key_text=pyfaustcov)  

- For matfaust the CI job is ``matfaust_test_code_coverage`` it uses MOCov to
  produce a coverage report as this [one](https://faustgrp.gitlabpages.inria.fr/-/faust/-/jobs/3186956/artifacts/coverage_html/index.html)
  (this link might be deleted in the future, in which case you might look at artifacts of recent ci
  [jobs](https://gitlab.inria.fr/faustgrp/faust/-/jobs)). The tests used for coverage calculation
  are pretty much the same as the ones in ``ctest_matlab`` and ``doctest_matfaust``.
  Look at the badge below to know what is the current coverage of matfaust tests. If it is not so good
  it is partially due to the fact that the code for mex float (single in matlab) is not auto-tested
  but there is no reason to think it would fail more often than the code for double which is tested
  (except for tests based on a strong accuracy of the results obtained for approximation algorithms).
  However it would be advisable to pursue adding tests for a larger coverage.

![matfaust test coverage](https://gitlab.inria.fr/faustgrp/faust/badges/hakim_branch/coverage.svg?job=matfaust_test_code_coverage&key_text=matfaustcov&key_width=90)

- C++ test coverage: the calculation of the coverage for C++ tests is
  integrated directly in the ``ctest`` CI job
 (the badge below indicates the coverage rate).

- About GPU code: GPU code is tested mainly in C++ tests but ``pyfaust.tests`` submodule
  provides the same tests on GPU and CPU (except for the CPU only API). matfaust is
  not so well tested on GPU code, another reason for its poor coverage rate.  Not all
  VMs are capable to run GPU tests (cf. [README](README.md) for requirements), the gitlab-runners
  must be tagged 'cuda' if they are capable of running the GPU code (this condition
  is verified dynamically in the concerned CI jobs).

- Gitlab badges: Gitlab provides a badge feature for indicating clearly the test coverage rate. The
  pyfaust/matfaust/C++ core coverage badges are displayed above and below.

![C++ test coverage](https://gitlab.inria.fr/faustgrp/faust/badges/hakim_branch/coverage.svg?job=ctest&key_text=c++cov)  

### Package tests and automatic releases

Several CI jobs (defined in ``.gitlab-ci.yml``) are responsible to build and
generate different packages.

#### Jobs for revision packages

They are triggered on each push of commits. The goal of the packages is to work
on the project without the need to produce a new release. They are named with the
SHA256 prefix of the corresponding git-commit. The secondary interest is to
avoid to find out too late (at release time) that one package doesn't build.

- Jobs for system packages which embed both pyfaust and matfaust.
    * ``pkg_linux``: generates Linux packages \*.rpm (Fedora, Redhat, Centos, etc.) and
      \*.deb Linux packages (Debian, Ubuntu, etc.).
    * ``pkg_win``: generates \*.exe Windows 10 package.
    * ``pkg_macos``: generates \*.pkg Mac OS X package.

- Jobs for wheel PIP packages for pyfaust only (one per system type):
    * ``pkg_linux_purepy_rev``,
    * ``pkg_win_purepy_rev``,
    * ``pkg_macos_purepy_rev``.

#### Jobs for release packages

They are triggered on each git-tag push. The goal is to produce a release of
all packages available of the tagged version. The same kind of packages as the
ones generated on revision/commit are generated with additionally anaconda packages.

The jobs are listed below. See description of equivalent jobs for revision packages.

- Jobs for system packages which embed both pyfaust and matfaust
    * ``pkg_linux_release``
    * ``pkg_linux_release_static``: does the same as the previous one except that
      matio library is built statically in the wrappers (no need to install matio
      separately on the system)
    * ``pkg_macos_release``
    * ``pkg_win_release``

- Jobs for wheel PIP packages for pyfaust only (one per system type):
    * ``pkg_linux_purepy_release``
    * ``pkg_macos_purepy_release``
    * ``pkg_win_purepy_release``

- Jobs for anaconda packages: ``conda_linux_pub``, ``conda_macosx_pub``,
    ``conda_win_pub`` that produce anaconda packages based on wheel packages.

- Jobs for extra feature packages.

    * Another version of Python is used to build the wrappers in these jobs:
      ``pkg_linux_purepy_release_extra_pyver``, ``pkg_win_purepy_release_extra_pyver``,
      ``pkg_macos_purepy_release_extra_pyver``. We maintain two versions of Python.

    * ``pkg_macos_purepy_release_torch_linked``,
      ``pkg_linux_purepy_release_torch_linked``,
      ``pkg_linux_purepy_release_openblaso`` build packages that provide features
      associated with tierce scientific libraries.

In order to avoid a non-workable package release there are several CI
jobs to tests the package after building and before uploading to package
 repositories. There are all ``test_`` prefixed and their names (suffixes)
 directly refer to previous jobs described.

The packages are then automatically uploaded on the Gitlab project package registry
(CI job ``upload_pkgs_to_registry``), on PyPI for wheel package (CI job ``pypi_pub``)
and on Anaconda (``conda_*_pub``).

### Documentation CI job

On release (git-tag push) documentations are automatically generated and put online
on Gitlab pages. The job responsible for that task is the classic ``pages`` job located
as usual in ``.gitlab-ci.yml`` root file.
