======= decay-book =======

Resources for the book *Finite Difference Computing with Exponential
Decay Models* by H. P. Langtangen, published by Springer, 2016. This
is an open access book where the author has all rights and publishes
various "ebook formats":
"http://hplgit.github.io/decay-book/doc/web/index.html" in the present
repository.

===== Organization of directories =====

The book is written in "DocOnce": "https://github.com/hplgit/doconce" using
the directory organization described in the repo
"setup4book-doconce":
"http://hplgit.github.io/setup4book-doconce/doc/web/index.html".

===== Nicknames for chapters =====

Directories in this repo apply a nickname as a short form of a chapter:

 * `alg`: Algorithms and implementations
 * `analysis`: Analysis
 * `genz`: Generalizations
 * `models`: Models
 * `softeng`: Software engineering

===== Key directories =====

 * The book: `doc/.src/book`
 * Individual chapters: `doc/.src./chapters` (text in `*.do.txt` files)
 * Source code for chapter `X`: `doc/.src/chapters/X/src-X` (`X` is chapter nickname)
 * Figures for chapter `X`: `doc/.src/chapters/X/fig-X`
 * Slides for chapter `X`: `doc/.src/chapters/slides-X/` (slides in `*.do.txt` files)

===== Published material =====

 * PDF and HTML books: `doc/pub` (gh-pages branch)
 * Source code (examples): `src`

===== How to modify slides =====

 o Go to the directory for the chapter: `doc/.src/chapters/X`
 o Go to the `slides-X` subdirectory
 o Modify the `*.do.txt` files as desired
 o Compile both the textbook material *and* the slides with `bash make.sh`
   in the parent directory `doc/.src/chapters/X`
