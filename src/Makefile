# Edit the following variables as needed
HDF_INSTALL = /usr/lib/x86_64-linux-gnu/hdf5/serial
EXTLIB = -L$(HDF_INSTALL)/lib
CC          = clang++
CFLAGS      = -O3 -std=c++14 -fopenmp=libomp
GL_FLAG		= -DUSE_GL=true
LIB         = -ldl -lsz -lz -lm -lyaml-cpp
GL_LIB		= -lGL -lX11 -lGLU -lGLEW -lglfw

INCLUDE   = -I$(HDF_INSTALL)/include
LIBSHDF   = $(EXTLIB) $(HDF_INSTALL)/lib/libhdf5_hl_cpp.a $(HDF_INSTALL)/lib/libhdf5_cpp.a $(HDF_INSTALL)/lib/libhdf5_hl.a $(HDF_INSTALL)/lib/libhdf5.a

protons:
	$(CC) $(CFLAGS) $(GL_FLAG) -o $@ main.cpp \
	                        config/*.cpp \
	                        detectors/*.cpp \
	                        fields/*.cpp \
	                        integrators/*.cpp \
	                        particles/*.cpp \
	                        sources/*.cpp \
	                        util/*.cpp \
                            graphics/*.cpp \
                            $(INCLUDE) $(LIBSHDF) $(LIB) $(GL_LIB)
	cp $@ ../runtime/$@

protons-headless:
	$(CC) $(CFLAGS) -o $@ main.cpp \
	                        config/*.cpp \
	                        detectors/*.cpp \
	                        fields/*.cpp \
	                        integrators/*.cpp \
	                        particles/*.cpp \
	                        sources/*.cpp \
                            util/*.cpp \
                            $(INCLUDE) $(LIBSHDF) $(LIB)
	cp $@ ../runtime/$@


clean: 
	rm -f *.o protons protons-headless
 
.SUFFIXES:.o.c

