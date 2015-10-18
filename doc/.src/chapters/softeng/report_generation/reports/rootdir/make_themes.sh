#!/bin/sh
# Make all themes given on the command line (or if no themes are
# given, make all themes in _themes/)

if [ $# -gt 0 ]; then
    themes=$@
else
    themes="ADCtheme agni alabaster basicstrap bloodish cbc fenics fenics_minimal1 fenics_minimal2 jal pylons scipy_lectures slim-agogo uio vlinux-theme agogo basic bizstyle classic default epub haiku nature pyramid scrolls sphinxdoc traditional bootstrap cloud redcloud  impressjs sphinx_rtd_theme"
fi

for theme in $themes; do
    echo "building theme $theme..."
    doconce subst -m "^html_theme =.*$" "html_theme = '$theme'" conf.py
    make clean
    make html
    mv -f _build/html sphinx-$theme
done
echo
echo "Here are the built themes:"
ls -d sphinx-*
echo "for i in sphinx-*; do google-chrome $i/index.html; done"

