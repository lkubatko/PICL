# Detect OS
ifeq ($(OS), Windows_NT)
    PLATFORM := windows
else
    UNAME := $(shell uname -s)
    ifeq ($(UNAME), Darwin)
        PLATFORM := mac
    else
        PLATFORM := linux
    endif
endif

# Output binary
TARGET := picl

# -ffp-contract=off disables fused multiply-add contraction. Without it, GCC
# is free to fuse a*b+c into a single FMA instruction whenever it helps, and
# whether it does so can depend on optimization level / surrounding code, not
# just on the source text. For most code that's an unnoticeable last-bit
# rounding difference; for this program's branch-length hill-climbing and
# annealing searches (bl_uphill_*/bl_anneal_* in anneal.c etc.), an
# accept/reject decision can sit close enough to a tie that a one-ULP change
# flips it, sending the search down a different (still statistically valid,
# but no longer reproducible-with-the-same-seed) trajectory. This flag keeps
# results reproducible across optimization levels/machines for a given seed.
FPFLAGS := -ffp-contract=off

ifeq ($(PLATFORM), mac)

    # Search for Boost in common locations
    BOOST_CANDIDATES := \
        $(wildcard ./boost_*/boost/math/special_functions/bessel.hpp) \
        $(wildcard $(HOME)/boost_*/boost/math/special_functions/bessel.hpp) \
        /opt/homebrew/include/boost/math/special_functions/bessel.hpp \
        /usr/local/include/boost/math/special_functions/bessel.hpp \
        /opt/local/include/boost/math/special_functions/bessel.hpp

    BOOST_HEADER := $(firstword $(foreach p, $(BOOST_CANDIDATES), $(wildcard $(p))))

    ifeq ($(BOOST_HEADER),)
        $(error Could not find Boost. Download it from https://www.boost.org/users/download/ \
            and extract it in this directory or your home directory)
    endif

    BOOST_INC := $(subst /boost/math/special_functions/bessel.hpp,,$(BOOST_HEADER))
    $(info Found Boost at: $(BOOST_INC))

    all:
	gcc -O2 -g $(FPFLAGS) -c main.c -o main.o
	g++ -O2 -g $(FPFLAGS) -c bessel_wrapper_macos.cpp -std=c++17 -I$(BOOST_INC) -o bessel_wrapper.o
	g++ main.o bessel_wrapper.o -lm -o $(TARGET)

    # Debug build: -O0 (no optimization, easier to trust line numbers), full
    # symbols, and AddressSanitizer + UndefinedBehaviorSanitizer, which turn
    # memory bugs (out-of-bounds access, use-after-free, etc.) into an
    # immediate, precise crash report with a source line, instead of a bare
    # "Segmentation fault". Produces picl_debug (much slower to run -- only
    # use it to diagnose a crash, not for real analyses).
    debug:
	gcc -O0 -g $(FPFLAGS) -fsanitize=address,undefined -fno-omit-frame-pointer -c main.c -o main_debug.o
	g++ -O0 -g $(FPFLAGS) -fsanitize=address,undefined -c bessel_wrapper_macos.cpp -std=c++17 -I$(BOOST_INC) -o bessel_wrapper_debug.o
	g++ -fsanitize=address,undefined main_debug.o bessel_wrapper_debug.o -lm -o $(TARGET)_debug

    clean:
	rm -f main.o main_debug.o bessel_wrapper.o bessel_wrapper_debug.o $(TARGET) $(TARGET)_debug

else ifeq ($(PLATFORM), windows)

    # bessel_wrapper.cpp uses std::cyl_bessel_k, a C++17 standard-library
    # function. -std=c++17 must be passed explicitly when compiling it --
    # without it, GCC silently falls back to an older default C++ standard
    # (which one depends on the GCC version) where cyl_bessel_k isn't
    # declared at all, and the build fails with "'cyl_bessel_k' is not a
    # member of 'std'". Compiling main.c and bessel_wrapper.cpp as separate
    # steps (like the mac branch already does) avoids passing a C++-only
    # flag to the plain-C compile of main.c.
    all:
	gcc -O2 -g $(FPFLAGS) -c main.c -o main.o
	g++ -O2 -g $(FPFLAGS) -std=c++17 -c bessel_wrapper.cpp -o bessel_wrapper.o
	g++ main.o bessel_wrapper.o -lm -lstdc++ -o $(TARGET)

    clean:
	del /f main.o bessel_wrapper.o $(TARGET).exe 2>nul || true

else

    # See the note in the windows branch above: -std=c++17 is required for
    # std::cyl_bessel_k to be visible, and isn't implied by default on every
    # GCC version. This is also the fix needed on HPC clusters like OSC and
    # OSU ASC's Unity, where the "cyl_bessel_k is not a member of std"
    # compile error commonly shows up -- if it persists after this fix, the
    # loaded/default compiler is likely too old to support C++17's special
    # math functions at all (GCC < 7); load a newer one first, e.g.
    # `module load gnu/<version>` (OSC) or the cluster's Lmod GCC module
    # (Unity), then re-run make.
    all:
	gcc -O2 -g $(FPFLAGS) -c main.c -o main.o
	g++ -O2 -g $(FPFLAGS) -std=c++17 -c bessel_wrapper.cpp -o bessel_wrapper.o
	g++ main.o bessel_wrapper.o -lm -lstdc++ -o $(TARGET)

    debug:
	gcc -O0 -g $(FPFLAGS) -fsanitize=address,undefined -fno-omit-frame-pointer -c main.c -o main_debug.o
	g++ -O0 -g $(FPFLAGS) -std=c++17 -fsanitize=address,undefined -c bessel_wrapper.cpp -o bessel_wrapper_debug.o
	g++ -fsanitize=address,undefined main_debug.o bessel_wrapper_debug.o -lm -lstdc++ -o $(TARGET)_debug

    clean:
	rm -f main.o bessel_wrapper.o main_debug.o bessel_wrapper_debug.o $(TARGET) $(TARGET)_debug

endif

.PHONY: all clean debug
