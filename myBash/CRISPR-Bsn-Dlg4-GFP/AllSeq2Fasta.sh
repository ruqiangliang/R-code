#/bin/sh
echo "Usage: > all.seq"
echo "This bash script combine all the .seq into one fasta format file"

for f in *.seq;
do 
  echo ">${f%.seq}";
  cat $f;
done
