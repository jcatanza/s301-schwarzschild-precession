"""pytest configuration: make the project's flat top-level modules
(orbit.py, campaign.py, ...) importable from tests/ without installing a
package, and force a non-interactive matplotlib backend so importing
fit_orbit.py never tries to open a display."""

import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
