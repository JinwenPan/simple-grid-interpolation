CC = g++
CFLAGS  = -O3 -Wall -Winline -Wshadow -std=c++17
TARGET = intug

all: main.cpp 
	$(CC) $(CFLAGS) -g main.cpp -o $(TARGET)