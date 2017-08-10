### Copyright (c) 2016 Luigi Rocca <luigi.rocca (@T) ge.imati.cnr.it>
### Released under the MIT License.

###///   -- -- -- -- -- eggMake -- -- -- -- --   \\\###
###    Library - compilation of a single target.    ###
###\\\   -- -- -- -- -- ------- -- -- -- -- --   ///###


### Include default values
egglib_SELF_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
include $(egglib_SELF_DIR)defaults.mk

### Set all directories that may contain the source files.
VPATH = $(egg_SOURCE_PATH) $(egg_DIR_SOURCES)

### Get filenames from directories
egglib_FULL_SOURCES := $(egg_SOURCES)
egglib_FULL_SOURCES += $(foreach dir,$(egg_DIR_SOURCES),$(wildcard $(dir)/*.c $(dir)/*.cpp $(dir)/*.cc $(dir)/*.cxx $(dir)/*.c++ $(dir)/*.C))
egglib_SOURCES := $(notdir $(egglib_FULL_SOURCES))

### Non-file phony targets
.PHONY : all clean release debug trace static tracerelease tracedebug staticrelease staticdebug statictrace statictracerelease statictracedebug



### Initialize flags
egglib_FINAL_FLAGS = $(CPPFLAGS) $(egg_FLAGS)

### Use release flags for `make', debug flags for `make debug' and
### tracing flags for "make trace".
release: egglib_FINAL_FLAGS += $(egg_RELEASE_FLAGS)
debug: egglib_FINAL_FLAGS += $(egg_DEBUG_FLAGS)
trace: egglib_FINAL_FLAGS += $(egg_TRACE_FLAGS)
static: egglib_FINAL_FLAGS += $(egg_STATIC_FLAGS)

tracerelease: egglib_FINAL_FLAGS += $(egg_TRACE_FLAGS) $(egg_RELEASE_FLAGS)
tracedebug: egglib_FINAL_FLAGS += $(egg_TRACE_FLAGS) $(egg_DEBUG_FLAGS)

staticdebug: egglib_FINAL_FLAGS += $(egg_STATIC_FLAGS) $(egg_DEBUG_FLAGS)
staticrelease: egglib_FINAL_FLAGS += $(egg_STATIC_FLAGS) $(egg_RELEASE_FLAGS)
statictrace: egglib_FINAL_FLAGS += $(egg_STATIC_FLAGS) $(egg_TRACE_FLAGS)
statictracerelease: egglib_FINAL_FLAGS += $(egg_STATIC_FLAGS) $(egg_TRACE_FLAGS) $(egg_RELEASE_FLAGS)
statictracedebug: egglib_FINAL_FLAGS += $(egg_STATIC_FLAGS) $(egg_TRACE_FLAGS) $(egg_DEBUG_FLAGS)

### Add specific language flags to general ones.
egglib_FINAL_CFLAGS = $(CFLAGS) $(egg_CFLAGS) $(egglib_FINAL_FLAGS)
egglib_FINAL_CXXFLAGS = $(CXXFLAGS) $(egg_CXXFLAGS) $(egglib_FINAL_FLAGS)

### Initialize linker flags
egglib_FINAL_LDFLAGS = $(LDFLAGS)
egglib_FINAL_LDLIBS = $(LDLIBS)

### dynamic
all: egglib_FINAL_LDFLAGS += $(egg_LDFLAGS)
all: egglib_FINAL_LDFLAGS += $(egg_LDLIBS)
release: egglib_FINAL_LDFLAGS += $(egg_LDFLAGS)
release: egglib_FINAL_LDLIBS += $(egg_LDLIBS)
debug: egglib_FINAL_LDFLAGS += $(egg_LDFLAGS)
debug: egglib_FINAL_LDLIBS += $(egg_LDLIBS)
trace: egglib_FINAL_LDFLAGS += $(egg_LDFLAGS)
trace: egglib_FINAL_LDLIBS += $(egg_LDLIBS)
tracerelease: egglib_FINAL_LDFLAGS += $(egg_LDFLAGS)
tracerelease: egglib_FINAL_LDLIBS += $(egg_LDLIBS)
tracedebug: egglib_FINAL_LDFLAGS += $(egg_LDFLAGS)
tracedebug: egglib_FINAL_LDLIBS += $(egg_LDLIBS)

### static
static: egglib_FINAL_LDFLAGS += $(egg_STATIC_LDFLAGS)
static: egglib_FINAL_LDLIBS += $(egg_STATIC_LDLIBS)
staticrelease: egglib_FINAL_LDFLAGS += $(egg_STATIC_LDFLAGS)
staticrelease: egglib_FINAL_LDLIBS += $(egg_STATIC_LDLIBS)
staticdebug: egglib_FINAL_LDFLAGS += $(egg_STATIC_LDFLAGS)
staticdebug: egglib_FINAL_LDLIBS += $(egg_STATIC_LDLIBS)
statictrace: egglib_FINAL_LDFLAGS += $(egg_STATIC_LDFLAGS)
statictrace: egglib_FINAL_LDLIBS += $(egg_STATIC_LDLIBS)
statictracerelease: egglib_FINAL_LDFLAGS += $(egg_STATIC_LDFLAGS)
statictracerelease: egglib_FINAL_LDLIBS += $(egg_STATIC_LDLIBS)
statictracedebug: egglib_FINAL_LDFLAGS += $(egg_STATIC_LDFLAGS)
statictracedebug: egglib_FINAL_LDLIBS += $(egg_STATIC_LDLIBS)

### Use C++ as a linker if there are C++ sources.
ifeq ($(strip $(filter %.cpp %.c++ %.cxx %.cc %.C, $(egglib_SOURCES)) ), )
egglib_LD = $(CC) $(egglib_FINAL_CFLAGS)
else
egglib_LD = $(CXX) $(egglib_FINAL_CXXFLAGS)
endif


ifeq ($(egg_TARGET_TYPE),executable)
egglib_TARGET := $(egg_TARGET)
else ifeq ($(egg_TARGET_TYPE),library)
egglib_TARGET := $(egg_LIB_PREFIX)$(egg_TARGET).$(egg_STATIC_LIB_E) $(egg_LIB_PREFIX)$(egg_TARGET).$(egg_DYNAMIC_LIB_E)
else ifeq ($(egg_TARGET_TYPE),static-library)
egglib_TARGET := $(egg_LIB_PREFIX)$(egg_TARGET).$(egg_STATIC_LIB_E)
else ifeq ($(egg_TARGET_TYPE),dynamic-library)
egglib_TARGET := $(egg_LIB_PREFIX)$(egg_TARGET).$(egg_DYNAMIC_LIB_E)
else
$(error Wrong value for egg_TARGET_TYPE variable)
endif


### Standard goals supported by eggmake:
all : $(egglib_TARGET)
release : $(egglib_TARGET)
debug : $(egglib_TARGET)
trace : $(egglib_TARGET)
static : $(egglib_TARGET)
tracerelease : $(egglib_TARGET)
tracedebug : $(egglib_TARGET)
staticrelease : $(egglib_TARGET)
staticdebug : $(egglib_TARGET)
statictrace : $(egglib_TARGET)
statictracerelease : $(egglib_TARGET)
statictracedebug : $(egglib_TARGET)


### Add automatically generated sources to the project sources
egglib_SOURCES += $(egg_AUTO_SOURCES)

### Set the list of modules to be compiled.
#ALL_EXTENSIONS := $(addprefix ".", $(CXX_EXTENSIONS)) .c
#OBJECTS := $(foreach ext, $(ALL_EXTENSIONS),   )
egglib_OBJECTS :=
egglib_OBJECTS += $(patsubst %.cpp, %.opp, $(filter %.cpp, $(egglib_SOURCES)) )
egglib_OBJECTS += $(patsubst %.c++, %.o++, $(filter %.c++, $(egglib_SOURCES)) )
egglib_OBJECTS += $(patsubst %.cxx, %.oxx, $(filter %.cxx, $(egglib_SOURCES)) )
egglib_OBJECTS += $(patsubst %.cc, %.oo, $(filter %.cc, $(egglib_SOURCES)) )
egglib_OBJECTS += $(patsubst %.C, %.O, $(filter %.C, $(egglib_SOURCES)) )
egglib_OBJECTS += $(patsubst %.c, %.o, $(filter %.c, $(egglib_SOURCES)) )
egglib_DEPS := 
egglib_DEPS += $(patsubst %.cpp, %.dpp, $(filter %.cpp, $(egglib_SOURCES)) )
egglib_DEPS += $(patsubst %.c++, %.d++, $(filter %.c++, $(egglib_SOURCES)) )
egglib_DEPS += $(patsubst %.cxx, %.dxx, $(filter %.cxx, $(egglib_SOURCES)) )
egglib_DEPS += $(patsubst %.cc, %.dd, $(filter %.cc, $(egglib_SOURCES)) )
egglib_DEPS += $(patsubst %.C, %.D, $(filter %.C, $(egglib_SOURCES)) )
egglib_DEPS += $(patsubst %.c, %.d, $(filter %.c, $(egglib_SOURCES)) )



### Build directory management.

# if build dir is non-empty:
ifneq ($(strip $(egg_BUILD_DIR)),)
egglib_BUILD_DIR_SLASH := $(egg_BUILD_DIR)/
ifneq ($(MAKECMDGOALS),clean)
$(shell if [ ! -d $(egg_BUILD_DIR) ]; then mkdir $(egg_BUILD_DIR); fi)
endif
endif
# prepend its name to targets.
override egglib_BUILD_DEPS := $(addprefix $(egglib_BUILD_DIR_SLASH), $(egglib_DEPS))
override egglib_BUILD_OBJECTS := $(addprefix $(egglib_BUILD_DIR_SLASH), $(egglib_OBJECTS))



### Rebuild everything from scratch if current make command is
### different from the previous one.

ifeq ($(egg_CMDGOAL_FORCEBUILD), true)
#$(warning egg_CMDGOAL_FORCEBUILD)

ifeq ($(MAKECMDGOALS),)
egglib_NEWGOAL := all
else
egglib_NEWGOAL := $(MAKECMDGOALS)
endif
egglib_OLDGOAL := $(shell if [ -f $(egg_CMDGOAL_FORCEBUILD_FILE) ]; then cat $(egg_CMDGOAL_FORCEBUILD_FILE); else echo __empty__; fi)
#$(warning $(OLDGOAL))
#$(warning $(NEWGOAL))
egglib_dummy := $(shell echo $(egglib_NEWGOAL) > $(egg_CMDGOAL_FORCEBUILD_FILE))

ifneq ($(egglib_NEWGOAL),clean)
ifneq ($(egglib_OLDGOAL),clean)
ifneq ($(MAKECMDGOALS),$(egglib_OLDGOAL))
ifneq ($(egglib_NEWGOAL), $(egglib_OLDGOAL))

#$(warning DELETE)
#$(warning $(RM) $(egglib_BUILD_OBJECTS) $(egglib_BUILD_DEPS) $(egglib_TARGET))
egglib_dummy := $(shell $(RM) $(egglib_BUILD_OBJECTS) $(egglib_BUILD_DEPS) $(egglib_TARGET) $(egg_AUTO_SOURCES) $(egg_AUTO_HEADERS))

endif
endif
endif
endif

endif



### Remove all compilation files.
clean :
ifneq ($(egg_CMDGOAL_FORCEBUILD_FILE),)
	-$(RM) $(egg_CMDGOAL_FORCEBUILD_FILE)
endif
# Warning, the following clean command is repeated in new_target_forcerebuild section.
	-$(RM) $(strip $(egglib_BUILD_OBJECTS) $(egglib_BUILD_DEPS) $(egglib_TARGET) $(egg_AUTO_SOURCES) $(egg_AUTO_HEADERS))
ifneq ($(strip $(egg_BUILD_DIR)),)
	-$(RM) -r $(strip $(egg_BUILD_DIR))
endif



### Add the makefile itself as a build dependency.
ifeq ($(egg_MAKEFILE_FORCEBUILD),true)
$(egglib_BUILD_DEPS) $(egglib_BUILD_OBJECTS) : $(firstword $(MAKEFILE_LIST))
endif


### Include dependencies, but only if we're not cleaning.
ifneq ($(MAKECMDGOALS),clean)
-include $(egglib_BUILD_DEPS)
endif



# ### Program linking:
# $(egg_TARGET) : $(egglib_BUILD_OBJECTS)
# 	$(strip $(egglib_LD) $^ $(egg_FRAMEWORKS) -o $@ $(egglib_FINAL_LDFLAGS) $(egglib_FINAL_LDLIBS))


### Linking and compilation commands

### Linking executable
ifeq ($(egg_TARGET_TYPE),executable)
define egglib_link_com
$(egglib_LD) $^ $(egg_FRAMEWORKS) -o $@ $(egglib_FINAL_LDFLAGS) $(egglib_FINAL_LDLIBS)
endef
$(egg_TARGET) : $(egglib_BUILD_OBJECTS)
	$(strip $(egglib_link_com))
endif
### Dynamic library
ifneq (,$(filter $(egg_TARGET_TYPE),library dynamic-library))
define egglib_dynamic_lib_com
$(egglib_LD) $^ $(egg_FRAMEWORKS) -o $@ $(egglib_FINAL_LDFLAGS) $(egglib_FINAL_LDLIBS)
endef
$(egg_LIB_PREFIX)$(egg_TARGET).$(egg_DYNAMIC_LIB_E) : $(egglib_BUILD_OBJECTS)
	$(strip $(egglib_dynamic_lib_com))
endif
### Static library
ifneq (,$(filter $(egg_TARGET_TYPE),library static-library))
define egglib_static_lib_com
$(AR) $(egg_ARFLAGS) $@ $^ 
endef
$(egg_LIB_PREFIX)$(egg_TARGET).$(egg_STATIC_LIB_E) : $(egglib_BUILD_OBJECTS)
	$(strip $(egglib_static_lib_com))
endif


define egglib_cc_dep
$(strip $(CC) $(egglib_FINAL_CFLAGS) -c $< -o $@)
endef

define egglib_cc_com
$(strip $(CC) -MM -MP -MG -MF $@ -MT $(egglib_BUILD_DIR_SLASH)$*.o -MT $@ $(egglib_FINAL_CFLAGS) $<)
endef

define egglib_cxx_dep
$(strip $(CXX) -MM -MP -MG -MF $@ -MT $(egglib_BUILD_DIR_SLASH)$*.$1 -MT $@ $(egglib_FINAL_CXXFLAGS) $<)
endef

define egglib_cxx_com
$(strip $(CXX) $(egglib_FINAL_CXXFLAGS) -c $< -o $@)
endef

### Implicit rules for objects and dependencies, C language.
$(egglib_BUILD_DIR_SLASH)%.o: %.c
	$(egglib_cc_dep)

$(egglib_BUILD_DIR_SLASH)%.d: %.c
	$(egglib_cc_com)

### Implicit rules for objects and dependencies, C++ language.
$(egglib_BUILD_DIR_SLASH)%.opp: %.cpp
	$(egglib_cxx_com)
$(egglib_BUILD_DIR_SLASH)%.oxx: %.cxx
	$(egglib_cxx_com)
$(egglib_BUILD_DIR_SLASH)%.o++: %.c++
	$(egglib_cxx_com)
$(egglib_BUILD_DIR_SLASH)%.oo: %.cc
	$(egglib_cxx_com)
$(egglib_BUILD_DIR_SLASH)%.O: %.C
	$(egglib_cxx_com)

$(egglib_BUILD_DIR_SLASH)%.$.dpp: %.cpp
	$(call egglib_cxx_dep,opp)
$(egglib_BUILD_DIR_SLASH)%.$.dxx: %.cxx
	$(call egglib_cxx_dep,oxx)
$(egglib_BUILD_DIR_SLASH)%.$.d++: %.c++
	$(call egglib_cxx_dep,o++)
$(egglib_BUILD_DIR_SLASH)%.$.dd: %.cc
	$(call egglib_cxx_dep,oo)
$(egglib_BUILD_DIR_SLASH)%.$.D: %.C
	$(call egglib_cxx_dep,O)


# $(egglib_BUILD_DIR_SLASH)%.$.dpp: %.cpp
# 	$(strip $(CXX) -MM -MP -MG -MF $@ -MT $(egglib_BUILD_DIR_SLASH)$*.opp -MT $@ $(egglib_FINAL_CXXFLAGS) $<)
# $(egglib_BUILD_DIR_SLASH)%.$.dxx: %.cxx
# 	$(strip $(CXX) -MM -MP -MG -MF $@ -MT $(egglib_BUILD_DIR_SLASH)$*.oxx -MT $@ $(egglib_FINAL_CXXFLAGS) $<)
# $(egglib_BUILD_DIR_SLASH)%.$.d++: %.c++
# 	$(strip $(CXX) -MM -MP -MG -MF $@ -MT $(egglib_BUILD_DIR_SLASH)$*.o++ -MT $@ $(egglib_FINAL_CXXFLAGS) $<)
# $(egglib_BUILD_DIR_SLASH)%.$.dd: %.cc
# 	$(strip $(CXX) -MM -MP -MG -MF $@ -MT $(egglib_BUILD_DIR_SLASH)$*.oo -MT $@ $(egglib_FINAL_CXXFLAGS) $<)
# $(egglib_BUILD_DIR_SLASH)%.$.D: %.C
# 	$(strip $(CXX) -MM -MP -MG -MF $@ -MT $(egglib_BUILD_DIR_SLASH)$*.O -MT $@ $(egglib_FINAL_CXXFLAGS) $<)



# ### Implicit rules for objects and dependencies, C++ language.
# ### Foreach version, still broken
# $(BUILD_DIR_SLASH)%.$(OBJ_EXT): %.$(CPP_EXT)
# 	$(CXXC) $(CXXDEFINES) $(CXXFLAGS) $(CXXINCLUDES) -c $< -o $@  

# $(BUILD_DIR_SLASH)%.$(DEP_EXT): %.$(CPP_EXT)
# 	$(CXXC) -MM -MP -MG -MF $@ -MT $(BUILD_DIR_SLASH)$*.$(OBJ_EXT) -MT $@ $(CXXDEFINES) $(CXXFLAGS) $(CXXINCLUDES) $<
# #	$(CXXC) -MM $(CXXDEFINES) $(CXXFLAGS) $(CXXINCLUDES) $< | sed 's/\($(subst /,\/,$*)\)\.o[ :]*/$(BUILD_DIR_SLASH_SED)\1.$(OBJ_EXT) $(subst /,\/,$@) : /g' > $@;


# ### Old version inside a define, again not working.
# define CXXRULES

# OBJ_EXT = $(subst c,o,$(1)) 
# DEP_EXT = $(subst c,d,$(1)) 

# $(BUILD_DIR_SLASH)%.$(OBJ_EXT): %.$(1)
# 	$(CXXC) $(CXXDEFINES) $(CXXFLAGS) $(CXXINCLUDES) -c $< -o $@  

# $(BUILD_DIR_SLASH)%.$(DEP_EXT): %.$(1)
# 	$(CXXC) -MM -MP -MG -MF $@ -MT $(BUILD_DIR_SLASH)$*.$(OBJ_EXT) -MT $@

# endef

# $(warning $(CXX_EXTENSIONS))
# $(warning $(CXXRULES))
# $(warning $(foreach ext, $(CXX_EXTENSIONS), $(eval $(call CXXRULES, $(ext))) ) )


# ### debugging facilities
# comment out to debug the makefile
# all: myprint
# myprint:
# ### debug print
#	echo $(OBJECTS)
#	echo $(DEPS)
# $(warning $(egglib_OBJECTS))
# $(error error)

