LIB_DAI := /home/stream/lib/libDAI/libDAI-0.2.3

RM := rm -rf

OBJS := \
./src/main/cpp/CBP.o \
./src/main/cpp/FactorGraphProb.o \
./src/main/cpp/CFactorGraph.o \
./src/main/cpp/Compress.o \
./src/main/cpp/PositionCompress.o \
./src/main/cpp/LiftedGibbs.o \

CPP_DEPS := \
./src/main/cpp/CBP.d \
./src/main/cpp/FactorGraphProb.d \
./src/main/cpp/CFactorGraph.d \
./src/main/cpp/Compress.d \
./src/main/cpp/PositionCompress.d \
./src/main/cpp/LiftedGibbs.d \

WRAPPER_OBJS := \
./src/main/cpp/wrapperDAI.o \
./src/main/cpp/wrapperSTREAM.o \

WRAPPER_CPP_DEPS := \
./src/main/cpp/wrapperDAI.d \
./src/main/cpp/wrapperSTREAM.d \

# All Target
all: libSTREAM.so libSTREAMWrapper.so

# rule for Python wrapper object
src/main/cpp/wrapper%.o: src/main/cpp/wrapper%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DDEMO_DOCSTRING_SHOW_ALL=true -I"./include" -I$(LIB_DAI)"/include" -I/usr/include/python2.6 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

libSTREAMWrapper.so: $(WRAPPER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L$(LIB_DAI)/lib -L"." -Wl,-no-undefined -shared -o"libSTREAMWrapper.so" $(WRAPPER_OBJS) -lboost_python -lutil -lpthread -ldl -ldai -lpython2.6 -lSTREAM
	@echo 'Finished building target: $@'
	@echo ' '

# rule for library objects
src/main/cpp/%.o: src/main/cpp/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DNO_FANCY_COMP_ALG -I"./include" -I$(LIB_DAI)"/include" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

libSTREAM.so: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L$(LIB_DAI)"/lib" -Wl,-no-undefined -shared -o"libSTREAM.so" $(OBJS) -ldai
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS) $(CPP_DEPS) $(WRAPPER_OBJS) $(WRAPPER_CPP_DEPS) libSTREAM.so libSTREAMWrapper.so
	-@echo ' '


