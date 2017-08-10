### Copyright (c) 2016 Luigi Rocca <luigi.rocca (@T) ge.imati.cnr.it>
### Released under the MIT License.

###///  -- -- -- -- -- eggMake -- -- -- -- --  \\\###
###    Library - support for multiple targets.    ###
###\\\  -- -- -- -- -- ------- -- -- -- -- --  ///###



egglib_PHONY_TARGETS := clean release debug trace static tracerelease tracedebug staticrelease staticdebug statictrace statictracerelease statictracedebug
.PHONY : all $(egglib_PHONY_TARGETS) $(egg_TARGET_LIST)

all: $(egg_TARGET_LIST)


# define egglib_PASSRULE
# ifneq ($(filter $(1), $(MAKECMDGOALS)),)
# egglib_PASSTARGET := $(1)
# ifeq ($(words $(MAKECMDGOALS)),1)
# $(1) : $(egg_TARGET_LIST)
# endif
# endif
# endef

define egglib_PASSVAR =
$(if $(filter $(1), $(MAKECMDGOALS)),$(eval egglib_PASSTARGET := $(1)))
endef

define egglib_PASSRULE =
$(if $(filter $(1), $(MAKECMDGOALS)),$(if $(filter $(words $(MAKECMDGOALS)),1),$(eval $(1) : $(egg_TARGET_LIST))))
endef

$(foreach rule, $(egglib_PHONY_TARGETS), $(call egglib_PASSVAR, $(rule)))
$(foreach rule, $(egglib_PHONY_TARGETS), $(call egglib_PASSRULE, $(rule)))



# egg_NO_PRINT_DIR defaults to true
ifeq ($(origin egg_NO_PRINT_DIR), undefined)
egg_NO_PRINT_DIR := true
endif

# Don't print the directory change when calling make again,
# if egg_NO_PRINT_DIR is true.
ifeq ($(egg_NO_PRINT_DIR),true)
egglib_NO_PRINT_DIR_OPT=--no-print-directory
else
egglib_NO_PRINT_DIR_OPT=
endif

$(egg_TARGET_LIST):
	@$(MAKE) $(egglib_NO_PRINT_DIR_OPT) -f Makefile.$@ $(egglib_PASSTARGET)


$(egglib_PHONY_TARGETS):
	@echo nothing > /dev/null
