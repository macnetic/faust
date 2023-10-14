# How to contribute

[[_TOC_]]

Below we describe the processes to follow ideally for any contribution to the FAµST
project.

### Posting issues on Gitlab

Please specify these elements on any issue which might need it:

- The system you're using (Linux, Mac OS X, Windows) for FAµST, its version.
- The wrapper you're using (pyfaust and matfaust).
- The installer you used (wheel package -- pip, conda, system packages -- .exe,
  .rpm, .deb, .pkg).
- The version of FAµST: check ``pyfaust.__version__`` on a Python terminal for
  pyfaust or ``matfaust.version()`` on a Matlab terminal for matfaust.
- Explain clearly the problem you ran into or the feature you're interesting in.
  For any bug please give a snippet of code to reproduce easily.

### Code contribution

Any contribution to the development is of course welcome. It might be a bug
fix, a new test, a feature or an algorithm implementation.

In order to contribute there are two alternatives:

1. If you have an Inria account or are a member for the FAµST Gitlab project
   please fork the project (a button is provided on the Gitlab front-page of the
   project). Then do a Merge Request. If you are member of the project you might
   also directly push your branch (giving a meaningful name) and do your merge
   request for this branch.

2. If you don't have have any account on Inria's Gitlab you can git-clone the
  project, do any modifications you want to do. Then send us the patch obtained
  with the command ``git diff > my.patch``. Please send to all email addresses
  listed in [contact section of the official web site](https://faust.inria.fr/contact/).
  If you ask us we might fork the project for you and add your account as
  a project member (it needs first that an Inria sponsor set your account in the
  [external account portal](https://external-account.inria.fr)). Once your account
  is set you can proceed as in first alternative.

Before sending us your patches or merge requests, please ensure as much as
possible that the project is still able to build and that related tests pass
(but we'll help you in that process if you find it too complicated).
For any further information about building please consult the main
[README](README.md) (an easy Docker way is provided to build pyfaust wrappers)
and [tests README](README-tests-ci-release.md).

### Documentation / tutorials

You might also be interested in contributing with a tutorial. So far we have used Jupyter
Notebook for pyfaust tutorials and Matlab Live Scripts for matfaust. You can send
us your tutorial, if we find it useful we will add it on the [website](https://faust.inria.fr).
