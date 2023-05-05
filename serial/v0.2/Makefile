
OBJ_DIR = debug

SOURCE_DIR = src

third_party = ../../thirdParty

SOURCE_FILES = $(SOURCE_DIR)/bcond.c \
	$(SOURCE_DIR)/geocalc.c \
	$(SOURCE_DIR)/ioutils.c \
	$(SOURCE_DIR)/main.c \
	$(SOURCE_DIR)/material.c \
	$(SOURCE_DIR)/mesh.c \
	$(SOURCE_DIR)/octree.c \
	$(SOURCE_DIR)/param.c \
	$(SOURCE_DIR)/post.c \
	$(SOURCE_DIR)/parser.c \
	$(SOURCE_DIR)/gradient.c

SOURCE_OBJ = $(OBJ_DIR)/bcond.o \
	$(OBJ_DIR)/geocalc.o \
	$(OBJ_DIR)/ioutils.o \
	$(OBJ_DIR)/main.o \
	$(OBJ_DIR)/material.o \
	$(OBJ_DIR)/mesh.o \
	$(OBJ_DIR)/octree.o \
	$(OBJ_DIR)/param.o \
	$(OBJ_DIR)/post.o \
	$(OBJ_DIR)/parser.o \
	$(OBJ_DIR)/gradient.o

# LASPACK
LASPACK_DIR = $(third_party)/LASPack 
LASPACK_LIB = $(third_party)/LASPack/laspack/liblaspack.a

EXECUTABLE = ../../bin/OpenFVM-serial-v0.2

CC = gcc
CFLAGS = -I$(SOURCE_DIR) -I$(LASPACK_DIR) -O2

all: $(SOURCE_OBJ) 

	$(CC) $(CFLAGS) $(SOURCE_OBJ) $(LASPACK_LIB) -o $(EXECUTABLE) -lm

$(OBJ_DIR)/bcond.o: $(SOURCE_DIR)/bcond.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/bcond.c -o $(OBJ_DIR)/bcond.o

$(OBJ_DIR)/geocalc.o: $(SOURCE_DIR)/geocalc.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/geocalc.c -o $(OBJ_DIR)/geocalc.o

$(OBJ_DIR)/ioutils.o: $(SOURCE_DIR)/ioutils.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/ioutils.c -o $(OBJ_DIR)/ioutils.o

$(OBJ_DIR)/main.o: $(SOURCE_DIR)/main.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/main.c -o $(OBJ_DIR)/main.o

$(OBJ_DIR)/material.o: $(SOURCE_DIR)/material.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/material.c -o $(OBJ_DIR)/material.o

$(OBJ_DIR)/mesh.o: $(SOURCE_DIR)/mesh.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/mesh.c -o $(OBJ_DIR)/mesh.o

$(OBJ_DIR)/octree.o: $(SOURCE_DIR)/octree.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/octree.c -o $(OBJ_DIR)/octree.o

$(OBJ_DIR)/param.o: $(SOURCE_DIR)/param.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/param.c -o $(OBJ_DIR)/param.o

$(OBJ_DIR)/post.o: $(SOURCE_DIR)/post.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/post.c -o $(OBJ_DIR)/post.o

$(OBJ_DIR)/parser.o: $(SOURCE_DIR)/parser.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/parser.c -o $(OBJ_DIR)/parser.o

$(OBJ_DIR)/gradient.o: $(SOURCE_DIR)/gradient.c
	$(CC) $(CFLAGS) -c $(SOURCE_DIR)/gradient.c -o $(OBJ_DIR)/gradient.o

clean:
	rm $(SOURCE_OBJ) $(LASPACK_OBJ)

