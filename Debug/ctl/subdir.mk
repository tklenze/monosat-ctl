################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ctl/CTLFormula.cpp 

OBJS += \
./ctl/CTLFormula.o 

CPP_DEPS += \
./ctl/CTLFormula.d 


# Each subdirectory must supply rules for building sources it contributes
ctl/%.o: ../ctl/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -I.././ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11  -Wno-unused-variable -Wno-unused-but-set-variable   -Wno-sign-compare -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


