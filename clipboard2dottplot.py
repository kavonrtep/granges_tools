#!/usr/bin/env python
"""
takes content of cliboard, save it to temporary fasta file and run dotter program on it
temporary fasta is deleted after dotter is closed
rrequires dotter and xsel programs to be installed
"""

import subprocess
import tempfile

# run xsel command and save output to variable
xsel = subprocess.Popen(['xsel', '-b'], stdout=subprocess.PIPE)
# save to fasta file
fasta = tempfile.NamedTemporaryFile(suffix='.fasta')
# include header
fasta.write(b'>clipboard\n')
fasta.write(xsel.stdout.read())
fasta.flush()

# run dotter
subprocess.call(['dotter', fasta.name, fasta.name])

# delete fasta file
fasta.close()

