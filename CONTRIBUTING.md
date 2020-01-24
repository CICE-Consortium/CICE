# Contributing to CICE

First, thanks a lot for considering to contribute to CICE !

## Getting help
If you still need help getting started using the model afer reviewing the [model documentation][doc-running], the first point of contact with the CICE Consortium is the [Consortium Community Forum][forum]. 
This forum is monitored by Consortium members and also opened to the whole community.
**Please do not use our issue tracker for general support questions.**

[doc-running]: https://cice-consortium-cice.readthedocs.io/en/master/user_guide/ug_running.html
[forum]: https://xenforo.cgd.ucar.edu/cesm/forums/cice-consortium.146/

## Reporting issues
We use GitHub issues mainly for Consortium-confirmed bugs in the model and discussion around Consortium-led future directions/development.
If you believe you have found a bug in CICE, please first [open a new thread][new-thread] on the [Forum][forum].
If you are not sure you're running into a bug, open a thread on the [Forum][forum].  
Be sure to document:
1. The steps you took
2. The expected behavior of the model
3. Its actual behavior
4. On which [supported machine](/configuration/scripts/machines) you are running CICE

If you are runnning CICE on your own machine, and have a porting-related problem, we suggest first seeking help from your local IT help desk as they might be more helpful than us for porting issues !
If you still need help, you can [open a new thread][new-thread] on the Forum.
Be sure to document:
1. Your operating system and version
2. The `env.<machine>_<compiler>` and `Macros.<machine>_<compiler>` you are using
Also, make sure to use approppriate [BB code syntax][bb-code] to format your post, especially if you include [preformatted text][bb-code-code].


For feature/enhancement proposals, we suggest first opening a [new thread][new-thread] on the [Forum][discussion-board-wiki].

[new-thread]: https://xenforo.cgd.ucar.edu/cesm/forums/cice-consortium.146/post-thread
[discussion-board-wiki]: https://github.com/CICE-Consortium/About-Us/wiki/Contacting-the-Consortium#discussion-board
[bb-code]: https://xenforo.cgd.ucar.edu/cesm/help/bb-codes/
[bb-code-code]: https://xenforo.cgd.ucar.edu/cesm/help/bb-codes/#code

## Contributing to development
New CICE features should first be discussed with the Consortium members in the [Forum][forum]. 

A quick start overview of how to get the code, compile the model and run a case is available in the [README][quick-start].

[quick-start]: README.md#getting-started

### Coding guidelines and developer documentation
Please review and follow our documented [Software Development Practices][dev-practices] when contributing code. 
To get started coding, you might want to read the [Developer Guide][developer-guide] section of the model documentation.

[dev-practices]: https://github.com/CICE-Consortium/About-Us/wiki/Software-Development-Practices
[developer-guide]: https://cice-consortium-cice.readthedocs.io/en/master/developer_guide/index.html

### Version control workflow
CICE uses a forking and Pull Request (PR) workflow, so start by [forking][cice-fork] the repo to your GitHub account. 
A detailed overview of our Git workflow, as well as a Git primer for less experienced users is found on the [About-Us wiki][git-workflow].
We recommend you review this page appropriately.  

We use a feature branch workflow, so please keep your pull requests self-contained, i.e. your pull request branch should be separate from `master` and should not contain unrelated changes.  

Note: the column portion of the model, [Icepack][icepack], is included in CICE as a [Git submodule][wiki-submodule].
As submodules can be considered a more advanced feature of Git that even experienced Git users might never have encountered, we recommend getting confortable with the concept.
Some reading materials, in increasing level of complexity:
- The [Submodules section][git-tower] of the Git Tower tutorial (introductory)
- The [Submodules chapter][git-book-submodule] of the [Pro Git][pro-git] book (in depth)
- The official Git documentation: [gitsubmodules][gitsubmodules], [git-submodule][git-submodule], [gitmodules][gitmodules] (reference)

[cice-fork]: https://github.com/CICE-Consortium/CICE/fork
[git-workflow]: https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance
[icepack]: https://github.com/CICE-Consortium/Icepack
[wiki-submodule]: https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance#submodules
[git-book-submodule]: https://git-scm.com/book/en/v2/Git-Tools-Submodules
[pro-git]: https://git-scm.com/book/en/v2
[git-tower]: https://www.git-tower.com/learn/git/ebook/en/command-line/advanced-topics/submodules
[gitsubmodules]: https://git-scm.com/docs/gitsubmodules
[git-submodule]: https://git-scm.com/docs/git-submodule
[gitmodules]: https://git-scm.com/docs/gitmodules

### Testing
Any contribution to CICE must be properly tested to ensure that any new feature does not considerably change the model outputs if the feature is turned off.

The required level of testing is described briefly on the [About-Us wiki][req-testing], and more in depth in the [Testing CICE][doc-testing] section of the model documentation.
A detailed walk-through is available at [End-To-End Testing Procedure][end-to-end].  

PR's introducing new feature should have corresponding tests added to the [base_suite](configuration/scripts/tests/base_suite.ts).

[req-testing]: https://github.com/CICE-Consortium/About-Us/wiki/Software-Development-Practices#required-testing
[doc-testing]: https://cice-consortium-cice.readthedocs.io/en/master/user_guide/ug_testing.html
[end-to-end]: https://cice-consortium-cice.readthedocs.io/en/master/user_guide/ug_testing.html#end-to-end-testing-procedure

### Documentation
The source of the CICE model documentation is in the [`doc/source`](doc/source) folder of the repository, and is written in [reStructuredText][rst] using [Sphinx][sphinx].
It is automatically deployed from updates to the `master` branch to [Read the Docs][rtd].  

New features should be documented, and we suggest that any pull request introducing new code features also update the model documentation accordingly.
The [About-Us wiki][wiki-doc] has information about:
- [Frequently Asked Questions][faq] about the documentation
- [reStructuredText syntax][wiki-rst]
- [Building the documentation locally][wiki-local]

Opening a pull request will trigger a build of the documentation at `https://external-builds.readthedocs.io/html/cice-consortium-cice/<pr-number>/index.html` in order for reviewers to look at any documentation changes.

[rst]: https://docutils.sourceforge.io/rst.html
[sphinx]: http://www.sphinx-doc.org/en/master/
[rtd]: https://cice-consortium-cice.readthedocs.io/
[wiki-doc]: https://github.com/CICE-Consortium/About-Us/wiki/Documentation-Workflow-Guidance
[faq]: https://github.com/CICE-Consortium/About-Us/wiki/Documentation-Workflow-Guidance#faqs
[wiki-rst]: https://github.com/CICE-Consortium/About-Us/wiki/Documentation-Workflow-Guidance#editing-rst-files
[wiki-local]: https://github.com/CICE-Consortium/About-Us/wiki/Documentation-Workflow-Guidance#using-sphinx
[wiki-pr]: https://github.com/CICE-Consortium/About-Us/wiki/Documentation-Workflow-Guidance#push-changes-back-to-the-repository

### Code review and merging process
Once opened, your pull request will be reviewed by members of the Consortium. 
You might be asked to polish your PR, add tests or add documentation before it is integrated in the `master` branch.
Please respond to the given feedback and ask any question you might have.
Once the reviewers are satisfied with the pull request, it will be [squashed-merged][squash] into the `master` branch.
Any forward development at this point should happen in a new branch created from from the updated `master`.

[squash]: https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-request-merges#squash-and-merge-your-pull-request-commits