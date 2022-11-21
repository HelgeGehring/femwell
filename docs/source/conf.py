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