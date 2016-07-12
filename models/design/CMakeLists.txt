cmake_minimum_required(VERSION 2.8)
project(design)

execute_process(COMMAND python ../script.py)
install(FILES ./build/trento-config.txt DESTINATION trento)
install(FILES ./build/vhlle-config.txt DESTINATION vhlle)
install(FILES ./build/escape-config.txt DESTINATION escape)
install(FILES ./build/stat-config.txt DESTINATION collector)


