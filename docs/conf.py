# Licensed under a 3-clause BSD style license - see licenses/LICENSE.rst
#
# Sphinx configuration for the StarKit documentation. The shared Astropy
# settings come from sphinx-astropy; install them with
# pip install -e '.[docs]'.

import datetime

from sphinx_astropy.conf.v1 import *  # noqa: F401,F403

import starkit

# -- Project information ------------------------------------------------------

project = 'starkit'
author = 'Wolfgang Kerzendorf, Tuan Do'
copyright = '{0}, {1}'.format(datetime.datetime.now().year, author)

# The short X.Y version and the full version.
version = starkit.__version__.split('-', 1)[0]
release = starkit.__version__

# -- General configuration ----------------------------------------------------

exclude_patterns += ['_templates', '**.ipynb_checkpoints']

extensions += [
    'nbsphinx',
    'sphinx.ext.mathjax',
]

nbsphinx_prolog = """
The notebook is available here:
https://github.com/starkit/starkit/tree/master/docs/{{ env.doc2path(env.docname, base=None) }}

----
"""

nbsphinx_execute = 'never'

# -- Options for HTML output --------------------------------------------------

html_theme_options = {
    'logotext1': 'Star',  # white,  semi-bold
    'logotext2': 'Kit',  # orange, light
    'logotext3': ':documentation'   # white,  light
    }

html_title = '{0} v{1}'.format(project, release)

# Output file base name for HTML help builder.
htmlhelp_basename = project + 'doc'

# -- Options for LaTeX and manual page output ---------------------------------

latex_documents = [('index', project + '.tex', project + ' Documentation',
                    author, 'manual')]

man_pages = [('index', project.lower(), project + ' Documentation',
              [author], 1)]

# -- Resolving issue number to links in changelog -----------------------------
github_issues_url = 'https://github.com/starkit/starkit/issues/'
