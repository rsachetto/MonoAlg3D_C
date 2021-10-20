

all: release

release:
	@./build.sh

debug:
	@./build.sh -d

simulator:
	@./build.sh simulator

simulator_debug:
	@./build.sh -d simulator

clean:
	@./build.sh -d clean; ./build.sh clean
