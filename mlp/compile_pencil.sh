OPENCV_DIR="/home/b/src-not-saved/opencv-2.4.5-build"

cd ..
make
cd mlp

gcc -std=c99 -c mlp_impl_glu.pencil.c -o mlp_impl_glu.pencil.o -lm -I../opencl -I../core -I/opt/local/include -I/usr/local/cuda/include/ -lstdc++ -lm  -L$OPENCV_DIR/lib/  -lboost_serialization -lOpenCL  -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video   -I$OPENCV_DIR/include/ -Werror=int-to-pointer-cast -I../opencl -I../core

gcc -std=c99 -c mlp_impl_kernels.pencil.c -o mlp_impl_kernels.pencil.o -lm -I../opencl -I../core -I/opt/local/include -I/usr/local/cuda/include/ -lstdc++ -lm  -L$OPENCV_DIR/lib/  -lboost_serialization -lOpenCL  -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video   -I$OPENCV_DIR/include/  -Werror=int-to-pointer-cast -I../opencl -I../core

g++ mlp_impl_glu.pencil.o mlp_impl_kernels.pencil.o -o mlp_impl_glu.pencil -lm test_mlp.o allocator.o ../opencl/errors.o -std=c++0x -I../opencl -I../core -I/opt/local/include -I/usr/local/cuda/include/ -lstdc++ -lm  -L$OPENCV_DIR/lib/  -lboost_serialization -lOpenCL  -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video   -I$OPENCV_DIR/include/   -Werror=int-to-pointer-cast -I../opencl -I../core

ulimit -s unlimited
