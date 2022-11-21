import os
import re
from datetime import date

release = re.sub('^v', '', os.popen('git describe').read().strip())
version = release

project = 'Femwell'
project_copyright = f'{date.today().year}, Femwell developers'

extensions = [
    'matplotlib.sphinxext.plot_directive'
]

plot_rcparams = {
    'savefig.bbox': 'tight'
}

html_theme = "sphinx_book_theme"

html_theme_options = {
    "path_to_docs": "docs/source",
    "repository_url": "https://github.com/HelgeGehring/femwell",
    "repository_branch": "main",
    "use_edit_page_button": True,
    "use_issues_button": True,
    "use_repository_button": True,
    "use_download_button": True,
}