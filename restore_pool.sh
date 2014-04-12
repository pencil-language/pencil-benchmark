#/bin/bash

cd build/pool

rm -rf *.jpg
rm -rf *.xml

echo
echo "Restoring *.jpg and *.xml files in build/pool."

for f in `ls *.7z`; do
	base_name=`basename $f .7z`
	git checkout $base_name &> /dev/null
done

cd -
