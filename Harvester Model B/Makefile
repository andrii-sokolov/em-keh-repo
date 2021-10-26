CC = g++ 
CFLAGS = -Wall -O2
TARGET = emkeh
SOURCE = src/main.cpp src/optimization_CD.cpp src/optimization_HB.cpp src/optimization_interp.cpp src/optimization_MSQR.cpp src/optimization_out.cpp src/optimization_RK.cpp
LDFLAGS = -lboost_system -lboost_iostreams -lboost_filesystem
MAGNET_TARGET = mag_field
MAGNET_SOURCE = magnet_src/main.cpp magnet_src/field.cpp magnet_src/progressbar.cpp magnet_src/harmonic_ballance.cpp 
MAGNET_LDFLAGS = -lpthread


all: $(SOURCE)
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS); 
	cp Data/RMS/* ./;
	./$(TARGET);
	
magnet: $(MAGNET_SOURCE)
	$(CC) $(CFLAGS) $(MAGNET_SOURCE) -o $(MAGNET_TARGET) $(MAGNET_LDFLAGS);
	
clear: 
	rm *.txt;
	rm $(TARGET);
	rm *.dat;
	rm $(MAGNET_TARGET);
