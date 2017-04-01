SRC_DIR=output-front
TAR_DIR=output-front-jpg

mkdir -p $TAR_DIR
cd $SRC_DIR
for filename in *.tga; do
    convert $filename ../$TAR_DIR/${filename%.*}.jpg
done
