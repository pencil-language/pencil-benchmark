OPENCV_DIR="/home/b/src-not-saved/opencv-2.4.5-build"

cd ..
make
cd mlp

gcc -std=c99 -c mlp_impl.pencil.c -o mlp_impl.pencil.o -lm -I../opencl -I../core -I/opt/local/include -I/usr/local/cuda/include/ -lstdc++ -lm  -L$OPENCV_DIR/lib/  -lboost_serialization -lOpenCL  -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video   -I$OPENCV_DIR/include/ -Wimplicit -Werror=implicit-function-declaration -Werror=int-to-pointer-cast -I../opencl -I../core

g++ mlp_impl.pencil.o -o mlp_impl.pencil -lm test_mlp.o allocator.o ../opencl/errors.o -std=c++0x -I../opencl -I../core -I/opt/local/include -I/usr/local/cuda/include/ -lstdc++ -lm  -L$OPENCV_DIR/lib/  -lboost_serialization -lOpenCL  -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video   -I$OPENCV_DIR/include/  -Wimplicit -Werror=implicit-function-declaration -Werror=int-to-pointer-cast -I../opencl -I../core
