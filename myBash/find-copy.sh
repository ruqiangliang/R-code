#This script is used for search some kind of files
# then cp to destination folder with original folder structure
find . -name '*.bs' -exec cp --parents \{\} /home/ruqiang/Desktop/myBash \;

