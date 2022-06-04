.. $Id: intro.rst 4471 2013-11-21 19:44:22Z renard $

.. include:: ../replaces.txt

.. highlight:: none

.. _dp-contribute:

How to contribute / Git repository on Savannah
==============================================

.. |saweb| replace:: Savannah
.. _saweb: https://savannah.gnu.org

.. |sawebg| replace:: Getfem on Savannah
.. _sawebg: https://savannah.nongnu.org/projects/getfem

.. |linktask| replace:: here
.. _linktask: https://savannah.nongnu.org/task/?group=getfem
			
.. |sawebgsrc| replace:: Getfem sources on Savannah
.. _sawebgsrc: http://git.savannah.nongnu.org/gitweb/?p=getfem.git;a=tree

.. |tfweb| replace:: transifex
.. _tfweb: https://www.transifex.com

.. |tfwebteam| replace:: Getfem translation team on Transifex
.. _tfwebteam: https://www.transifex.com/tkoyama010/getfem-doc/dashboard

.. |cfvlang| replace:: Currently supported languages by Sphinx are
.. _cfvlang: https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-language

.. |sphintl| replace:: Sphinx Internationalization
.. _sphintl: http://www.sphinx-doc.org/en/master/intl.html

.. |readthedocs| replace:: Read the Docs
.. _readthedocs: https://getfem.readthedocs.io

|gf| is an  open source finite element library based on a collaborative development. If you intend to make some contributions, you can ask for membership of the project there. Contributions of all kinds are welcome: documentation, bug reports, constructive comments, changes suggestions, bug fix, new models, etc ...

Contributors are of course required to be careful that their changes do not affect the proper functioning of the library and that these changes follow a principle of backward compatibility.

See |linktask|_ for a list of task and discussions about |gf| development.

**IMPORTANT** : a contributor implicitly accepts that his/her contribution will be distributed under the LGPL licence of |gf|.

The main repository of |gf| is on Savannah, the software forge of the Free Software Foundation (see |saweb|_). The page of the project on Savannah is |sawebg|_. See also |sawebgsrc|_.

How to get the sources
----------------------

.. |sagit| replace:: git on Savannah
.. _sagit: http://savannah.gnu.org/maintenance/UsingGit/

If you just want the sources and do not intend to make some contributions, you can just use the command ::

  git clone https://git.savannah.nongnu.org/git/getfem.git

If you intend to make some contributions, the first step is to ask for the inclusion in the |gf| project (for this you have to create a Savannah account). You have also to register a ssh key (see |sagit|_) and then use the command ::

  git clone ssh://savannah-login@git.sv.gnu.org:/srv/git/getfem.git

How to contribute
-----------------

Before modifying any file, you have to create a *development branch* because it is *not allowed to make a modification directly in the master branch*. It is recommended that the branch name is of the type `devel-name-subject` where name is your name or login and subject the main subject of the changes. For instance, if you chose `devel-me-rewrite-fem-kernel` as the branch name, the creation of the branch reads ::

  git branch devel-me-rewrite-fem-kernel
  git checkout devel-me-rewrite-fem-kernel

The first command create the branch and the second one position you on your branch. After that you are nearly ready to makes some modifications. You can specify your contact name and e-mail with the following commands in order to label your changes ::

  git config --global user.name "Your Name Comes Here"
  git config --global user.email you@yourdomain.example.com


Specific branch for doc improvements and typo-fixes
---------------------------------------------------

If you want to contribute to the documentation only, it is not necessary to build a specific branch. You can just checkout to the ``fixmisspell`` branch which has been created for this purpose with ::

  git checkout fixmisspell


Locally commit your changes
---------------------------

Once you made some modifications of a file or you added a new file, say `src/toto.cc`, the local commit is done with the commands::

  git add src/toto.cc
  git commit -m "Your extensive commit message here"

At this stage the commit is done on your local repository but not in the Savannah one.

Push you changes in the Savannah repository
-------------------------------------------

You can now transfer your modifications to the Savannah repository with ::

  git push origin devel-me-rewrite-fem-kernel

where of course *devel-me-rewrite-fem-kernel* is still the name of your branch. At this stage your modifications are registered in the branch *devel-me-rewrite-fem-kernel* of Savannah repository.
Your role stops here, since you are not allowed to modify the master branch of |gf|.


Ask for an admin to merge your modifications to the master branch of |gf|
-------------------------------------------------------------------------

Once you validated your modifications with sufficient tests, you can ask an admin of |gf| to merge your modifications. For this, contact one of them directly, or send an e-mail to *getfem-commits@nongnu.org* with the message : "please merge branch devel-me-rewrite-fem-kernel" with eventually a short description of the modifications. IMPORTANT : by default, your branch will be deleted after the merge, unless you express the need to keep it.


Merge modifications done by other contributors
----------------------------------------------

You can run a ::

  git pull origin master
  git merge master

in order to integrate the modifications which has been validated and integrated to the master branch. This is recommended to run this command before any request for integration of a modification in the master branch.


Some useful git commands
------------------------
::

  git status  : status of your repository / branch

  git log --follow "filepath"   : Show all the commits modifying the specified file (and follow the eventual change of name of the file).

  gitk --follow filename : same as previous but with a graphical interface


Contributing to document translation
------------------------------------

The recommended way for new contributors to translate document is to join |tfwebteam|_ . For contribution, please make account in |tfweb|_ and click request language and fill form . After translation, pull translated po file from site by using transifex-client. You need api token which you can get in transifex site. ::

  cd doc/sphinx
  tx pull -l <lang>

Set code for your native language to <lang> (see |cfvlang|_ ).

.. warning::

  **DO NOT** tx push to transifex. It will have some trouble. You can upload file one by one in team page.

After pulling translated po files, set <lang> to LANGUAGE in `doc/sphinx/Makefile.am` . ::

  LANGUAGE      = <lang>
  SPHINXOPTS    = -D language=$(LANGUAGE)

Then, you can run a following commands in order to make html localization document. ::

  cd doc/sphinx
  make html

If you want to make pdf file in your language, you can run a ::

  make latex
  cd build/latex
  make all-pdf-<lang>

See details in |sphintl|_ .

You can see translated document at |readthedocs|_ by switch language.
