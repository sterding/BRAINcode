#!/bin/sh

ls -lh $1 | awk '{gsub(/\/accepted_hits.bam/,"",$9); print $9}' | awk '{gsub(/\/[\/a-zA-Z0-9\_\/]*\//,"", $1); print $1}' > t1.txt
ls -lh $1 | awk '{print $9}' > t2.txt
wc  t1.txt > tlen.txt
read lines words characters filename < tlen.txt

seq $lines | awk '{print "sample"$1}' > t3.txt

paste t1.txt t2.txt t3.txt > t4.txt

echo -e Sample'\t'IDBam'\t'FileNotes > thead.txt

cat thead.txt t4.txt > samplelist_file.txt

rm t1.txt t2.txt t3.txt t4.txt t5.txt tlen.txt thead.txt
