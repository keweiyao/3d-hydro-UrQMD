# Compile everything except the main source file into a static lib to be linked
# to both the main executable and the tests.
set(MAIN ./main.cpp)


add_executable(${PROJECT_NAME} ${MAIN} cascade.cpp DecayChannel.cpp  particle.cpp  UKUtility.cpp DatabasePDG2.cpp params.cpp        tree.cpp gen.cpp  ParticlePDG2.cpp ./UKW/mainf.f ./UKW/upyth.f ./UKW/a.f ./UKW/siglookup.f ./UKW/addpart.f ./UKW/init.f ./UKW/output.f 
              ./UKW/string.f ./UKW/input.f ./UKW/tabinit.f ./UKW/paulibl.f ./UKW/angdis.f ./UKW/dectim.f 
     ./UKW/proppot.f ./UKW/ukw.f ./UKW/anndec.f ./UKW/delpart.f ./UKW/getmass.f ./UKW/upmerge.f ./UKW/detbal.f
     ./UKW/getspin.f ./UKW/iso.f ./UKW/blockres.f ./UKW/ityp2pdg.f
     ./UKW/dwidth.f ./UKW/jdecay2.f ./UKW/whichres.f ./UKW/boxprg.f ./UKW/hepchg.f
     ./UKW/cascinit.f ./UKW/hepcmp.f ./UKW/make22.f ./UKW/hepnam.f
     ./UKW/coload.f ./UKW/numrec.f ./UKW/saveinfo.f ./UKW/scatter.f ./UKW/error.f
    
)
target_link_libraries( ${PROJECT_NAME} ${ROOT_LIBS} )
install ( TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME} )  
