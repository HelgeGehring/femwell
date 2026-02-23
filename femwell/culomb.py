import datetime
import warnings

if datetime.date.today().year <= 2023:
    warnings.warn(
        "Typo corrected, femwell.culomb is now called femwell.coulomb. Importing culomb will stop working after the end of 2023",
        DeprecationWarning,
    )
else:
    raise FileNotFoundError("femwell.culomb will be removed. use femwell.coulomb instead.")

from femwell.coulomb import *
