NAME = main
SRCDIR = ./source

# RAYLIB = ./libs/raylib-5.0_win64_mingw-w64
CFILES := $(wildcard $(SRCDIR)/*.c)

LINUX_CC = gcc
LINUX_OBJDIR = ./objects/linux
LINUX_OFILES = $(patsubst $(SRCDIR)/%.c, $(LINUX_OBJDIR)/%.o, $(CFILES))
LINUX_RAYFLAGS = -lraylib -lGL -lm -lpthread -ldl -lrt -lX11

# WIN_CC = x86_64-w64-mingw32-gcc
# WIN_OBJDIR = ./objects/windows
# WIN_OFILES = $(patsubst $(SRCDIR)/%.c, $(WIN_OBJDIR)/%.o, $(CFILES))

OFLAGS = -O0
CFLAGS = -Wall -g -I$(SRCDIR) $(OFLAGS) 

# x86_64-w64-mingw32-gcc source/main.c -o main.exe -I ./libs/raylib-5.0_win64_mingw-w64/include/ -L ./libs/raylib-5.0_win64_mingw-w64/lib/ -l :libraylib.a -mwindows -lwinmm
# WIN_LINK_FLAGS = -L$(RAYLIB)/lib -l:libraylib.a -mwindows -lwinmm -lm
# WIN_INCL_FLAGS = -I$(RAYLIB)/include/

all: linux
	@echo "\033[1;33m ~> Compiling for both linux and windows"

linux: $(NAME)

$(NAME): $(LINUX_OFILES)
	$(LINUX_CC) -o $(NAME) $^ $(LINUX_RAYFLAGS)

$(LINUX_OBJDIR)/%.o:$(SRCDIR)/%.c
	$(LINUX_CC) $(CFLAGS) -c -o $@ $<

# win: $(NAME).exe

# $(NAME).exe: $(WIN_OFILES)
# 	$(WIN_CC) -o $(NAME).exe $^ $(WIN_INCL_FLAGS) $(WIN_LINK_FLAGS) 

# $(WIN_OBJDIR)/%.o:$(SRCDIR)/%.c
# 	$(WIN_CC) $(CFLAGS) -o $@ -c $(WIN_INCL_FLAGS) $<

clean:
	@echo "\033[1;31m ~> Cleaning objects files and executables\033[1;0m"
	-rm $(LINUX_OBJDIR)/*.o main

git:
	make clean
	@echo "\033[1;33m ~> Adding to git\033[1;0m"
	git add .
	@echo "\033[1;32m ~> Making commit to git\033[1;0m"
	git commit
	@echo "\033[1;34m ~> Pushing to github\033[1;0m"
	git push
