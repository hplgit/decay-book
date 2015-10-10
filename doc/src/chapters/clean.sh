#!/bin/sh -x
# To be run from a chapter directory
if [ -f papers.pub ]; then
    echo "Must be run from chapter directory, not this parent directory!"
    exit 1
fi

doconce clean

rm -rf .ptex2tex.cfg* newcommands*.tex automake_sphinx.py _*.html *.html *~ slides-*/*~ .*.exerinfo tmp* _doconce_deb* *-beamer*.pdf *-4print.pdf *-4screen.pdf *-plain.* \#* libpeerconnection.log *.pdf deck.js reveal.js latex_figs html_images *.tex _minted-*

if [ -d src-$1 ]; then
  cd src-$1
  if [ -f clean.sh ]; then
    sh clean.sh
  fi
  rm -f *.pyc
fi
cd ..
