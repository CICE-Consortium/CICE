:tocdepth: 3 

.. _doc:

Documentation System
====================

With CICE development, corresponding updates or modification to the CICE
documentation are required. Whenever you modify the model you should update
documentation. CICE uses `readthedocs.org <readthedocs.org>`_ to create 
online HTML and PDF documentation.

FAQs
----

1) What is reStructuredText (RST)?

   The CICE and Icepack documentation is written using reStructuredText (RST) markup language. 
   ReStructuredText is a markup language, like HTML, markdown, or LaTeX. 
   `readthedocs.org <readthedocs.org>`_ is a tool for publishing RST documents in other formats 
   such as HTML and PDF. Additional information about using RST and `readthedocs.org <readthedocs.org>`_ 
   are found in the sections below.

2) What is expected of *me* when changing the documentation?

   Updated static PDF documentation will be generated for each new CICE code release. However, 
   the online "master" version of HTML or PDF documentation is considered a living document and 
   will be updated regularly with regular code development workflow. 

   We expect that if you need to add or modify documentation that you will be able to modify the 
   RST source files and generate HTML in order to review the HTML documentation. We 
   will review the RST and HTML during a Pull Request to verify it is working properly and is consistent 
   with the rest of the CICE-Consortium documentation format. Then we will trigger a new documentation build
   on `CICE's documentation page <https://readthedocs.org/projects/cice-consortium-cice/>`_ when 
   the Pull Request is successful. The new documentation build will create the HTML and PDF versions of
   CICE's documentation along with your updates. 

   In particular, it is important that you test out tables, equations, section references, figures, and/or citations
   in your contributed documentation as these can be particularly fiddly to get right.


3) Where are the documentation files kept?

   The RST source files for generating HTML and PDF are stored in the master branch of the repository under /doc/source/. 

   The HTML and PDF versions of the documentation are available at `CICE's
   documentation page <https://readthedocs.org/projects/cice-consortium-cice/>`_
   HTML documentation for the current "master" branch as well as static documentation for releases of CICE 
   will be available on the `Versions page <https://readthedocs.org/projects/cice-consortium-cice/versions/>`_ 
   while corresponding PDF documentation is available on the `Downloads page
   <https://readthedocs.org/projects/cice-consortium-cice/downloads/>`_. The CICE-Consortium team will trigger
   builds of both HTML and PDF documentation with each pull request. 

.. _moddocs:

Steps for Modifying Documentation
---------------------------------

Setting up readthedocs.org
~~~~~~~~~~~~~~~~~~~~~~~~~~

The CICE-Consortium recommends that developers use `readthedocs.org <readthedocs.org>`_ to generate and test
their contributions to the CICE documentation. This tool does not require external libraries to be built
on each developer's personal computer and is free and easy to use. You can follow the steps below and also
reference the `Getting Started <https://docs.readthedocs.io/en/latest/getting_started.html>`_ guide available from `readthedocs.org <readthedocs.org>`_. 

1. Sign up for a free account at `readthedocs.org <readthedocs.org>`_

   Select a username and password. These do not have to match your GitHub username and password, but having
   the same username can be simpler if the user choses to do this. Below, 
   USERNAME is a placeholder - you would need to replace this with your personal username. 

2. Connect your GitHub account

   Click on your username in the upper right hand corner and select 'Settings' and then select 'Connected
   Services' to connect your GitHub account. This process will ask you to authorize a connection to
   readthedocs.org that allows for reading of information about and cloning of your repositories.

3. Import your projects

   Click on your username in the upper right hand corner and select 'My Projects'. Then click the 'Import
   a Project' green button. This will generate a list of repositories you are able to import. To add a
   repository click the + button. Once added and imported to readthedocs, the icon will change to an 
   upward pointing arrow. 

4. Modify the project settings

   Click on the project you are editing then click on the 'Admin' tab on the far right. The CICE-Consortium
   has found the following settings to be important for proper building.

   Under 'Settings' modify and save the following:
      - Name: USERNAME CICE    (this is the local name of the repository on readthedocs.org)
      - Repository URL: https://github.com/USERNAME/CICE.git
      - Repository type: Git
      - Description: (Add anything that would be useful to you to describe your fork.)
      - Documentation type: Sphinx html
      - Language: English
      - Programming Language: Only Words

   Under 'Advanced Settings' modify the following:
      - Install Project: Unchecked box
      - Requirements file: doc/requirements.txt  (*VERY IMPORTANT, see below*)
      - Default branch: readthedocs  (whatever branch you are working on for development. If not set, this will default to master.)
      - Default version: latest (what your documentation build will be called)
      - Enable PDF build: Checked box
      - Enable EPUB build: Checked box
      - Privacy Level: Public  (this is useful to keep public if you want to point to the tested documentation as part of a Pull Request)
      - Python Interpreter: Python 2.x


Model sandbox and documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Follow the general `CICE-Consortium Git Workflow and Developer's guide <https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance>`_
to clone the repository and create your personal fork for model modifications. Whenever you modify the model 
you should update documentation. You can update the documentation on the same branch of your fork on which 
you test code, or you can create a separate branch called 'readthedocs' to test only the RST and HTML documentation.

There are some important files you will need in order to correctly build the documentation. These should all be included automatically when you fork from the CICE-Consortium repositories:

   - /doc/requirements.txt : This file is necessary to get the references and citations working properly by using sphinxcontrib-bibtex. This file should *not* need to be modified by developers generally.
   - /doc/source/conf.py : Basic documentation information for the Consortium including template, etc. This file should *not* need to be modified by developers generally.
   - /doc/source/zreferences.rst : required for the references to link properly. This file should *not* need to be modified by developers generally. 
   - /doc/source/master_list.bib : the master list of references cited in the documentation. This file *may need* to be modified by developers with documentation updates. This file is currently ordered sequentially from oldest to newest and alphabetically within a given year. To add references for your documentation, edit the master_list.bib file using the Articles and/or Books entries as examples for your addition(s). Please follow the format for ordering the date/alphabetization as well as including a URL with the document's DOI.


Editing RST files
~~~~~~~~~~~~~~~~~~

Open the RST file using a text editor and make the changes necessary. Note that from the User's Guide documentation (see link above) there is a hyperlink called "Show Source" on the left hand column that will show you the RST source code for the HTML you are viewing. This is a good way to see the syntax for tables, equations, linking references, labeling tables or figures, and correctly identifying documentation sections or subsections.

Here are some resources for using RST files:

* `RST Primer1 <http://www.sphinx-doc.org/en/stable/rest.html>`_

* `RST Primer2 <http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_

* `RST Syntax <https://wiki.typo3.org/ReST_Syntax>`_

* `RST tables <http://www.sphinx-doc.org/en/stable/rest.html#tables>`_ - Note that tables can be tricky in Sphinx and we prefer using `comma separated tables <http://docutils.sourceforge.net/docs/ref/rst/directives.html#csv-table>`_ whenever possible.

Building documentation
~~~~~~~~~~~~~~~~~~~~~~

Once you've committed and pushed changes to the documentation `*.rst` files on your personal development fork. 
Go to your readthedocs.org site and then select your project "Overview". Whenever you commit to your fork
the documents will automatically build. There is also an option to "Build a Version". Choose "latest" 
and then click the green "Build version" button. 

You will automatically be taken to the "Builds" page with a list of recent documentation builds. 
The documentation build you just started will be listed as "Triggered" and then "Building". 
If the build is successful the status will change to "Passed"; if the build is not successful 
then the status will change to "Failed". You can click on the "Passed" or "Failed" text to get 
information about the build and what might be problematic. The time of the build attempt is also 
listed with the most recent build appearing at the top of the list.

To see the HTML you just successfully built, go to "Overview" and click on "latest" under versions. To see the PDF you just successfully built, go to "Downloads" and click on "latest PDF". 


Push changes back to the repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you're happy with the documentation you've generated, follow the standard CICE-Consortium 
`Git Workflow and Developer's guide <https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance>`_
to do a Pull Request and make sure to note in the Pull Request Template that documentation has also 
been updated. We will test the HTML and PDF as part of the Pull Request before it is merged to the repository. 
It can be particularly helpful if you include the link to your successfully built documentation that is 
part of the Pull Request, and in order to do this you must ensure that your settings in readthedocs.org 
are set to "Public".


Other Tips and Tricks
---------------------

Converting LaTeX to RST
~~~~~~~~~~~~~~~~~~~~~~~

If you start from a LaTeX (``*.tex``) document you will need to convert this to the RST format that Sphinx 
requires. A handy tool to do this is `Pandoc <http://pandoc.org/getting-started.html>`_, which you 
can install quickly and run from the command line.

Once Pandoc is installed, the basic command line syntax to convert a file is ::

     $ pandoc NAMEIN.tex -f latex -t rst -s -ou NAMEOUT.rst

The NAMEOUT.rst file can be directly edited for Sphinx. Pandoc does a beautiful job of converting the text, 
equations, and many tables. However, equation numbering, section linking, references, figures, and some 
tables required more hands on care to be sure they render correctly. 

Pandoc requires that the ``*.tex`` files be in utf-8 encoding. To easily do this open the ``*.tex``
document in Emacs then do ``ctrl-x ctrl-m f`` and you will be prompted to enter encoding type. Just
type in ``utf-8`` and hit enter. Then save with ``ctrl-x ctrl-s`` . You are done and the document can be
converted with Pandoc.

Using Sphinx
~~~~~~~~~~~~

We recommend that you use `readthedocs.org <readthedocs.org>`_ to test documentation
(see :ref:`moddocs`). However, it is also possible to use Sphinx to build and test documentation. 
If you choose to follow this workflow, below are some tips for using Sphinx. 

Installing Sphinx
`````````````````

Sphinx must be installed once on each platform. See `Sphinx <http://www.sphinx-doc.org/en/stable/>`_ or 
`Installing Sphinx <http://www.sphinx-doc.org/en/stable/install.html>`_ for details. Below are the
commands for installing Sphinx on a mac laptop at the command line. 
Other platforms may require other steps. ::

   $ sudo pip install --ignore-installed sphinx
   $ sudo pip install --ignore-installed sphinxcontrib-bibtex

The CICE Consortium has used the following software to get successful Sphinx HTML builds, including linked
references:

* python 2.7.11

* Sphinx (1.6.3)

* sphinx-rtd-theme (0.1.9)

* sphinxcontrib-bibtex (0.3.5)

* sphinxcontrib-websupport (1.0.1)

As mentioned above, you will need the conf.py, zreferences.rst, and master_list.bib files that are part of the 
master branch and automatically included in your checkout. To use linked references you will need to have the sphinxcontrib-bibtex package as well.

Building HTML
`````````````

Move into the /doc/ directory of your sandbox. Then execute the following command::

   $ make clean 

to get rid of old HTML files. Then execute::

   $ make html

to build HTML into /build/html/ directory. It will also give you errors if there is a problem with the build that will help you figure out how you need to modify your RST files for a successful HTML build. Finally ::

   $ open /build/html/FILE.html 

Open the HTML on your browser for testing.


Converting RST to PDF
`````````````````````

Generating a PDF is more complex and currently requires a two-step process. The generation will require 
recent versions of both LaTeX and Sphinx. From the /doc/ directory do the following::

     $ make latex
     $ cd build/latex
     $ make

Then search for the ``*.pdf`` document created.


