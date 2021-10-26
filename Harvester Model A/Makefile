CC = g++ 
CFLAGS = -Wall -O2
TARGET = emkeh
SOURCE = src/main.cpp src/rot_interp.cpp src/rot_minsqr.cpp src/rot_optimization.cpp src/rot_output.cpp src/rot_read_data.cpp src/rot_rk.cpp src/rot_semi_analytical.cpp src/rot.cpp src/transp_interp.cpp src/transp_minsqr.cpp src/transp_optimization.cpp src/transp_output.cpp src/transp_read_data.cpp src/transp_rk.cpp src/transp_semi_analytical.cpp src/transp.cpp
LDFLAGS = -lboost_system -lboost_iostreams -lboost_filesystem
MAGNET_TARGET = mag_field
MAGNET_SOURCE = magnet_src/main.cpp magnet_src/force.cpp
MAGNET_LDFLAGS = -lpthread


all: $(SOURCE)
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS); 
	cp Data/RMS/*.dat ./;
	./$(TARGET);
	
magnet: $(MAGNET_SOURCE)
	$(CC) $(CFLAGS) $(MAGNET_SOURCE) -o $(MAGNET_TARGET) $(MAGNET_LDFLAGS);
	
clear: 
	rm -f *.txt;
	rm -f $(TARGET);
	rm -f *.dat;
	rm -f $(MAGNET_TARGET);
