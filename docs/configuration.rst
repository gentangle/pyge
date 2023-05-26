Continuous Integration of the Project
======================================

This project builds on the Python Project template (`template`_) by Joao M. C. Teixeira. You can find more details on the project structure in its public repository. Below, I present the versioning scheme that follows the `Semantic
Versioning 2 <https://semver.org/>`_.

Every time a Pull Request is merged to `main` branch, the `deployment workflow
<https://github.com/joaomcteixeira/python-project-skeleton/blob/master/.github/workflows/version-bump-and-package.yml>`_
triggers. This action bumps the new version number according to the
merge commit message, creates a new GitHub tag for that commit, and
publishes in PyPI the new software version.

If the Pull Request merge commit starts with ``[MAJOR]``, a major version
increase takes place (attention to the rules of SV2!), if a PR merge commit
message starts with ``[FEATURE]`` it triggers a *minor* update. Finally, if the
commit message as not special tag, a *patch* update is triggered. Whether to
trigger a *major*, *minor*, or *patch* update concern mainly the main
repository maintainer.

This whole workflow can be deactivate if the commit to the ``main`` branch
starts with ``[SKIP]``.

In conclusion, every commit to ``main`` without the ``[SKIP]`` tag will be
followed by a version upgrade, new tag, new commit to ``main`` and consequent
release to PyPI. You have a visual representation of the commit workflow in the
`Network plot
<https://github.com/joaomcteixeira/python-project-skeleton/network>`_.

.. _template: https://github.com/joaomcteixeira/python-project-skeleton
