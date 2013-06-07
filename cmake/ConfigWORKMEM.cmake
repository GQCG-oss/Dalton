set(WORK_MEM_WORDS "64000000" CACHE STRING "DALTON WORK memory in words")
add_definitions(-DINSTALL_WRKMEM=${WORK_MEM_WORDS})

set(2EL_MEM_WORDS "1" CACHE STRING "DALTON static memory for storing 2-el integrals")
add_definitions(-DINSTALL_MMWORK=${2EL_MEM_WORDS})
