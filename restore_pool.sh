#/bin/bash

cd build/pool

rm -rf *.jpg
rm -rf *.xml

echo "Restoring *.jpg and *.xml files in build/pool."

for f in `ls *.7z`; do
	base_name=`basename $f .7z`
	echo $base_name
	git checkout $base_name
done

cd -
