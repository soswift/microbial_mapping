#!/bin/bash
for p in *.pdf
do 
   pdftoppm "$p" "${p%*.}" -png
done
