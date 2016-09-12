

ROOTLIBS = -L`root-config --libdir` -lMinuit `root-config --libs`
ROOTCFLAGS = `root-config --cflags`

CFLAGS = $(ROOTCFLAGS) -I.
LIBS = $(ROOTLIBS) -L. 

CC = g++ 

all: storehist_mee.exe fiteff_mee.exe storehist_fullmc.exe fiteff_mee_keys.exe

fiteff_mee_keys.exe: fiteff_mee_keys.o config.o
	$(CC) fiteff_mee_keys.o config.o -o fiteff_mee_keys.exe $(CFLAGS) $(LIBS)

fiteff_mee_keys.o: fiteff_mee_keys.cpp
	$(CC) -c fiteff_mee_keys.cpp -o fiteff_mee_keys.o $(CFLAGS)

storehist_fullmc.exe: storehist_fullmc.o config.o
	$(CC) storehist_fullmc.o config.o -o storehist_fullmc.exe  $(CFLAGS) $(LIBS) 

storehist_fullmc.o: storehist_fullmc.cpp
	$(CC) -c storehist_fullmc.cpp -o storehist_fullmc.o  $(CFLAGS)

fiteff_mee.exe: fiteff_mee.o config.o
	$(CC) fiteff_mee.o config.o -o fiteff_mee.exe $(CFLAGS) $(LIBS)

fiteff_mee.o: fiteff_mee.cpp
	$(CC) -c fiteff_mee.cpp -o fiteff_mee.o $(CFLAGS)

storehist_mee.exe: storehist_mee.o config.o
	$(CC) storehist_mee.o config.o -o storehist_mee.exe  $(CFLAGS) $(LIBS) 

config.o: config.cpp
	$(CC) -c config.cpp -o config.o 

storehist_mee.o: storehist_mee.cpp
	$(CC) -c storehist_mee.cpp -o storehist_mee.o  $(CFLAGS) 

clean: 
	\rm -f *.o
	\rm -f *.exe 
