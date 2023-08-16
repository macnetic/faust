# How to contribute

Below we describe the processes to follow for any contribution to the FAµST
project.

1. Posting issues on Gitlab
========================

Please specify these elements on any issue which might need it:

- The system you're using (Linux, Mac OS X, Windows) for FAµST, its version.
- The wrapper you're using (pyfaust and matfaust).
- The installer you used (wheel package -- pip, conda, system packages -- .exe,
  .rpm, .deb, .pkg).
- The version of FAµST: check ``pyfaust.__version__`` on a Python terminal for
  pyfaust or ``matfaust.version()`` on a Matlab terminal for matfaust.
- Explain clearly the problem you ran into or the feature you're intersting in
  and preferably give a snippet of code to reproduce easily.

2. Code contribution
====================

Any contribution to the development is of course welcome. It might be a bug
fix, a feature or an algorithm implementation.

In order to contribute:

- If you have an Inria account or are a member for the FAµST gitlab project
  please fork the project (a button is provided on the gitlab frontpage of the
  project). Then do a Merge Request. If you are member of the project you might
  also directly push your branch (giving a meaningful name) and do your merge
  request for this branch.

- If you don't have have any account on Inria's gitlab you can clone the
  project, do any modifications you want to do. Then send us the patch obtained
  with the command ``git diff > my.patch``. Please send to all email addresses
  listed in [contact section of the official web site](https://faust.inria.fr/contact/).

Before sending us your patches or merge requests, please ensure as much as
possible that the project is still able to build and that related tests pass
(but we'll help you in that process if you find it too complicated).
For any further information about building please consult the main
(README)[README.md] (an easy Docker way is provided to build pyfaust wrappers)
and (README-tests-ci-release.md) about tests.
