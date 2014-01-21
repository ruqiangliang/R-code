#/bin/sh
for f in *.xls;
do 
  ssconvert -S "$f" "${f%.xls}.csv";
done
