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
	gcc -g -c main.c -o main.o
	g++ -g -c bessel_wrapper_macos.cpp -std=c++17 -I$(BOOST_INC) -o bessel_wrapper.o
	g++ main.o bessel_wrapper.o -lm -o $(TARGET)

    clean:
	rm -f main.o bessel_wrapper.o $(TARGET)

else ifeq ($(PLATFORM), windows)

    all:
	gcc -g main.c bessel_wrapper.cpp -lm -lstdc++ -o $(TARGET)

    clean:
	del /f main.o bessel_wrapper.o $(TARGET).exe 2>nul || true

else

    all:
	gcc -g main.c bessel_wrapper.cpp -lm -lstdc++ -o $(TARGET)

    clean:
	rm -f main.o bessel_wrapper.o $(TARGET)

endif

.PHONY: all clean
