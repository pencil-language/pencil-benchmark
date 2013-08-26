OPENCV_DIR="/home/b/src-not-saved/opencv-2.4.5-build"

cd ..
make
cd mlp

gcc -std=c99 -c mlp_impl_glue.pencil.c -o mlp_impl_glue.pencil.o -lm -I../include -I../core -I/opt/local/include -I/usr/local/cuda/include/ -lstdc++ -lm  -L$OPENCV_DIR/lib/  -lboost_serialization -lOpenCL  -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video   -I$OPENCV_DIR/include/ -Werror=int-to-pointer-cast

gcc -std=c99 -c mlp_impl_kernels.pencil.c -o mlp_impl_kernels.pencil.o -lm -I../include -I../core -I/opt/local/include -I/usr/local/cuda/include/ -lstdc++ -lm  -L$OPENCV_DIR/lib/  -lboost_serialization -lOpenCL  -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video   -I$OPENCV_DIR/include/  -Werror=int-to-pointer-cast

g++ test_mlp.cpp serialization.cpp mlp_impl_glue.pencil.o mlp_impl_kernels.pencil.o -o ../mlp_impl_glue.pencil -lm -std=c++0x -I../include -I../core -I/opt/local/include -I/usr/local/cuda/include/ -lstdc++ -lm  -L$OPENCV_DIR/lib/  -lboost_serialization -lboost_system -lboost_filesystem -lOpenCL  -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video   -I$OPENCV_DIR/include/   -Werror=int-to-pointer-cast -L/usr/lib/

ulimit -s unlimited
