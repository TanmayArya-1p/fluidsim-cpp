build:
	cmake -S . -B build -DC‐MAKE_BUILD_TYPE=Release
	cd build && make

run : build
	./build/fluidsim

clean:
	rm -rf build/*

.PHONY: build clean run
