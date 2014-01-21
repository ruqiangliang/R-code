#/bin/sh

echo this is a bash script to split all .xls files in a folder into 
echo csv files, each sheet in one .xls file will be one csv file.
echo You need Gnumeric installed to do this.
for f in *.xls;
do 
  ssconvert -S "$f" "${f%.xls}.csv";
done
