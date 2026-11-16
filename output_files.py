# SPDX-License-Identifier: GPL-3.0-or-later
"""Écriture atomique des fichiers générés."""

import os
import tempfile
from contextlib import contextmanager


@contextmanager
def atomic_output(path):
    """Réserve un emplacement temporaire unique sur le même système de fichiers."""
    parent = os.path.dirname(os.path.abspath(path))
    with tempfile.TemporaryDirectory(dir=parent) as directory:
        temporary = os.path.join(directory, "output")
        yield temporary
        os.replace(temporary, path)
