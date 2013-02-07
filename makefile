

all: Makefile
	make -f Makefile
Makefile: 
	cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DITK_DIR=$(ITK_DIR)
clean:
	rm -rf CMakeCache.txt Makefile CMakeFiles/ ITKIOFactoryRegistration/ cmake_install.cmake  ITKIOFactoryRegistration/
tags:
	ctags -R --langmap=c++:+.txx --langmap=c++:+.cl $(ITK_SOURCE) .
