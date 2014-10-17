# I used this to cmpile theCython code. Locations of libraries etc should be changed.
# It's a huge pain so if you really need to recompile it... good luck!

g++ -c -O2 -fpic -std=c++0x corca.cpp
g++ -shared -o libcorca.so corca.o
mkdir build
mkdir build/temp.linux-x86_64-2.7
cython -a orca.pyx
g++ -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -O2 -fPIC -I. -I/project/soft/linux64/epd_free-7.3-2-rh5-x86_64/lib/python2.7/site-packages/numpy/core/include -I/project/soft/linux64/epd_free-7.3-2-rh5-x86_64/include/python2.7 -c orca.c -o build/temp.linux-x86_64-2.7/orca.o
g++ -pthread -shared -g -L. build/temp.linux-x86_64-2.7/orca.o -L/project/soft/linux64/epd_free-7.3-2-rh5-x86_64/lib -lcorca -lpython2.7 -o /project/home/ms7210/test_ABC/cython/orca/orca.so
