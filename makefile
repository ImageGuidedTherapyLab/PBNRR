

all: Makefile
	make -f Makefile
Makefile: 
	cmake -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.4 -DCMAKE_CXX_FLAGS_DEBUG='-g -O0' -DCMAKE_C_COMPILER=/usr/bin/gcc-4.4 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=ON -DITK_DIR=$(ITK_DIR)

clean:
	rm -rf CMakeCache.txt Makefile CMakeFiles/ ITKIOFactoryRegistration/ cmake_install.cmake  
tags:
	ctags -R --langmap=c++:+.txx --langmap=c++:+.cl $(ITK_SOURCE) .

.PHONY: tags
