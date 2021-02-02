#!/usr/bin/env bash

#(a) Change the current directory to a subdirectory called somedir
cd ./somedir

#(b) Print out to the terminal the contents of a file called sometext.txt.
cat ./sometext.txt

#(c) Print out to the terminal the last 5 lines of sometext.txt.
tail -n5 ./sometext.txt

#(d) Print out to the terminal the last 5 lines of each file that ends in the extension .txt
tail -n5 ./*.txt

#(e) Write a for loop which prints each integer from 0 to 6
for n in {0..6}; do echo $n; done
