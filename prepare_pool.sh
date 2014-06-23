cd build/pool
echo
echo "Removing *.jpg and *.xml files in build/pool"
rm -rf *.jpg
rm -rf *.xml
echo "Now extracting 1 image file and 1 xml file for tests"
git checkout M104_ngc4594_sombrero_galaxy_hi-res.jpg
cd -
