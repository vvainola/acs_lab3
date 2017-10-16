# --- MACROS

# define name of executable
EXEC = InfOli 

# define the C compiler to use
CC= gcc

# define any compile-time flags
CFLAGS= -O3
LDFLAGS = -I $(CUDA_INSTALL_PATH)/include -L $(CUDA_INSTALL_PATH)/lib64

# define any libraries to link into executable
ifeq ($(shell uname),Darwin)
	LDFLAGS = -framework OpenCL
else
	LDFLAGS += -lOpenCL -lm 
endif

# define C source files
SRCS= infoli.c ioFile.c init.c

# define C header files
HDRS= infoli.h kernel.h variables.h

# --- TARGETS
all: ${EXEC}

#Builds the program
${EXEC}: ${SRCS} ${HDRS}
	@echo #
	@echo "-- BUILDING PROGRAM --"
	${CC} ${SRCS} ${CFLAGS} ${LDFLAGS}  -o ${EXEC}

clean:
	@echo #
	@echo "-- CLEANING PROJECT FILES --"
	$(RM) *.o ${EXEC}

run: 
	./${EXEC}

RUNSERVER_REQ_RUN_SETTINGS=True
RUNSERVER_DEPS=all
-include $(ACS_SHARED_PATH)/runserver.mk
