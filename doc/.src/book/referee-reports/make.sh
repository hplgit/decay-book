#!/bin/sh
name=response
doconce format pdflatex $name --latex_code_style=vrb
pdflatex $name
