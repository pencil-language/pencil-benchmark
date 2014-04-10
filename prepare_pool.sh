cd build/pool
echo "Removing *.jpg and *.xml files in build/pool"
rm -rf *.jpg
rm -rf *.xml
echo "Now extracting 1 image file and 1 xml file for tests"
git checkout Kings_Cross_Western_Concourse_-_central_position_-_2012-05-02.75.jpg
cd -
