


INC_DIR = include
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

MKDIR = mkdir -p
BIN = $(BIN_DIR)/metcorr.exe 

SRC = $(wildcard $(SRC_DIR)/*.cc)
OBJ = $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

CC = g++

CPPFLAGS += -I$(INC_DIR)
CFLAGS +=  `root-config --cflags`
LDFLAGS += -Llib 
LDLIBS += `root-config --libs` -lMinuit2 
DEP += $(wildcard $(INC_DIR)/*.h)

.PHONY: all clean directories

all:  directories $(BIN)

directories: $(OBJ_DIR) $(BIN_DIR)

$(OBJ_DIR):
	$(MKDIR) $(OBJ_DIR)

$(BIN_DIR):
	$(MKDIR) $(BIN_DIR)

$(BIN): $(OBJ) 
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(DEP) 
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@ 



clean:
	$(RM) $(OBJ) ${BIN}
